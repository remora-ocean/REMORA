#ifndef REMORA_ICETHERMO_H_
#define REMORA_ICETHERMO_H_

#include <cmath>

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <algorithm>

#include <REMORA_SeaIce.H>

namespace REMORASeaIce {

/**
 * One column's worth of ice thermodynamics: a direct port of the body of
 * ROMS ice_thermo_tile (ROMS/Nonlinear/SeaIce/ice_mk.h), Mellor and Kantha
 * (1989) with Parkinson and Washington (1979).
 *
 * The routine is per-column on purpose. ROMS's kernel has no horizontal
 * stencil at all except the stress magnitude -- sustr(i,j) with sustr(i+1,j)
 * and svstr(i,j) with svstr(i,j+1) -- so that average is taken by the caller
 * and passed in as `utau_in`. Everything else here reads and writes a single
 * (i,j), which is why it can be an inline device function called from one
 * ParallelFor.
 *
 * Three pieces of ROMS's routine are deliberately absent, all of them dead in
 * the ICE_BULK_FLUXES configuration this port targets:
 *   - `d1`, `d2i`, `d3` (bulk transfer coefficients) are computed and never
 *     read; they belong to the path where the ice surface temperature is
 *     solved here rather than seeded by the bulk-flux routine.
 *   - `wsm` is set to zero and never read, and `t2` is computed twice and
 *     never read -- both are 2D scratch in ROMS, not outputs.
 *   - the first assignment to `Qi2` is overwritten before use.
 * Nothing else is omitted, reordered or simplified.
 */
struct IceColumn {
    // --- ocean, read ---
    amrex::Real temp_top;   //!< surface-level temperature [C]
    amrex::Real salt_top;   //!< surface-level salinity [psu]; the kernel clamps it to 0..40
    amrex::Real dztop;      //!< z_w(N) - z_r(N), surface half-cell thickness [m]

    // --- atmosphere, read ---
    amrex::Real Qnet_ai;    //!< net heat flux, ice to atmosphere [W/m2]
    amrex::Real Qnet_ao;    //!< net heat flux, open water to atmosphere [W/m2]
    amrex::Real snow;       //!< snowfall rate [m/s]
    amrex::Real rain;       //!< rainfall rate; ROMS treats rain < 0 as frozen
    amrex::Real utau_in;    //!< sqrt of the stress magnitude, averaged to rho

    // --- prognostic ice state, read and written ---
    amrex::Real Aice;       //!< concentration [fraction]
    amrex::Real Hice;       //!< thickness, cell average [m]
    amrex::Real Hsno;       //!< snow thickness [m]
    amrex::Real Hmel;       //!< melt-water thickness [m]
    amrex::Real Iage;       //!< age [days]
    amrex::Real Tice;       //!< interior temperature [C]

    // --- internal fields, read and written ---
    amrex::Real Isst;       //!< ice/snow surface temperature [C]
    amrex::Real T0mk;       //!< molecular sublayer temperature [C]
    amrex::Real S0mk;       //!< molecular sublayer salinity [psu]
    amrex::Real Qcon;       //!< heat conductivity, seeded by bulk_fluxes
    amrex::Real Qrhs;       //!< RHS heat flux, seeded by bulk_fluxes
    amrex::Real IOmf;       //!< ice-ocean mass flux [m/s]
    amrex::Real W_ai;       //!< melt/freeze rate, air/ice [m/s]
    amrex::Real W_ao;       //!< melt/freeze rate, air/ocean [m/s]
    amrex::Real W_fr;       //!< frazil accretion rate [m/s]; zero unless do_frazil
    amrex::Real W_io;       //!< melt/freeze rate, ice/ocean [m/s]
    amrex::Real W_ro;       //!< runoff rate [m/s]

    // --- ocean forcing, written ---
    amrex::Real stflx_temp; //!< surface heat flux [C m/s], ROMS sign convention
    amrex::Real stflx_salt; //!< surface salt flux, accumulated
};

/**
 * The one horizontal operation in ROMS's thermodynamics: the stress magnitude
 * at rho points, from ice_mk.h's
 *
 *     utau = SQRT(SQRT((0.5*(sustr(i,j)+sustr(i+1,j)))**2 +
 *                      (0.5*(svstr(i,j)+svstr(i,j+1)))**2))
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
ice_utau (const amrex::Real sustr_i, const amrex::Real sustr_ip1,
          const amrex::Real svstr_j, const amrex::Real svstr_jp1)
{
    const amrex::Real su = amrex::Real(0.5) * (sustr_i + sustr_ip1);
    const amrex::Real sv = amrex::Real(0.5) * (svstr_j + svstr_jp1);
    return std::sqrt(std::sqrt(su * su + sv * sv));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
ice_thermo_column (IceColumn& c,
                   const SeaIceParameters& p,
                   const amrex::Real dtice,
                   const amrex::Real rmask)
{
    using namespace Thermo;
    using amrex::Real;

    // ROMS fixes these as parameters regardless of rho0; p.use_rho0 routes
    // them to the model's own rho0 and Cp instead. See SeaIceParameters.
    const Real RhoSW  = p.RhoSW_eff();
    const Real Cp_w   = p.Cp_w_eff();
    const Real RhoCpr = p.RhoCpr_eff();

    const Real sec2day = Real(1.0) / Real(86400.0);

    // Extract ocean surface temperature and salinity.
    const Real temp_top = c.temp_top;
    const Real salt_top = std::min(std::max(Real(0.0), c.salt_top), Real(40.0));

    // Squared-root of surface wind stress magnitude.
    const Real utau = std::max(c.utau_in, Real(1.0e-4));

    // Snow thickness, ice thickness, brine fraction, thermal conductivity.
    const Real sice = std::min(Real(sice_ref), salt_top);
    const Real ice_thick  = Real(0.05) + c.Hice / (c.Aice + eps);
    const Real snow_thick = c.Hsno / (c.Aice + eps);
    Real brnfr = frln * sice / (c.Tice - eps);
    brnfr = std::min(brnfr, Real(0.2));
    brnfr = std::max(brnfr, Real(0.0));
    Real alph = AlphIc * std::max(Real(1.0) - Real(1.2) * brnfr, Real(0.25));
    {
        const Real cff = (c.Hice / Real(1.0)) * (c.Hice / Real(1.0));
        const Real corfac = Real(1.0) / (Real(0.5) * (Real(1.0) + std::exp(-cff)));
        alph *= corfac;
    }
    const Real Coa = Real(2.0) * alph * snow_thick / (AlphSn * ice_thick);

    // Solve for temperature at the top of the ice layer. Qcon and Qrhs arrive
    // seeded by bulk_fluxes with the linearisation of the air-ice flux.
    const Real b2d = Real(2.0) * alph / (ice_thick * (Real(1.0) + Coa));
    c.Qcon += b2d;
    if (c.Aice > p.min_ai) {
        c.Qrhs += b2d * (c.Tice + Real(273.15));
        c.Isst = (c.Qrhs / c.Qcon) - Real(273.15);
        // Bounded at zero Celsius for stability; it can occasionally take
        // Isst in the wrong direction.
        c.Isst = std::min(std::max(c.Isst, Real(-45.0)), Real(0.0));
    } else {
        c.Isst = temp_top;
    }

    // New interior ice temperature. The "cot" denominator is kept below zero
    // so that Tice can still move as it approaches zero Celsius.
    if (c.Aice > p.min_ai) {
        const Real tmin = std::min(c.Tice, frln * sice_ref);
        const Real cot  = -frln * sice * hfus / (tmin * tmin) + Cp_i;
        const Real cff1 = p.IceRho * cot * ice_thick * ice_thick;
        const Real cff2 = c.Isst - (Real(2.0) + Coa) * c.Tice;
        const Real cff3 = Real(1.0) + Coa;
        c.Tice += dtice * (Real(2.0) * alph / cff1 * (c.T0mk + cff2 / cff3));
        c.Tice = std::max(c.Tice, Real(-35.0));
        // Ice present: keep Tice below the reference freezing point, and below
        // the warmer of the surface and the water beneath.
        c.Tice = std::min(c.Tice, frln * sice_ref);
        c.Tice = std::min(c.Tice, std::max(c.Isst, c.T0mk));
    } else {
        c.Tice = temp_top;
    }

    // Associated heat fluxes.
    Real Qio = Real(0.0);
    if (c.Aice > p.min_ai) {
        const Real hicehinv = Real(1.0) / (Real(0.5) * ice_thick);
        Qio = alph * (c.T0mk - c.Tice) * hicehinv;
    }
    Real Qai = c.Qnet_ai;

    // Open water: set ice fluxes to zero.
    if (c.Aice <= p.min_ai) {
        c.Isst = c.T0mk;
        c.Tice = Real(-2.0);
        c.Isst *= rmask;
        c.Tice *= rmask;
        Qai  = Real(0.0);
        Qio  = Real(0.0);
        c.Hsno = Real(0.0);
        c.Hmel = Real(0.0);
    }

    // Surface water accumulation: ice melting.
    const Real ws = std::max(c.snow, Real(0.0));

    const Real tfrz = frln * sice;
    c.W_ai = Real(0.0);
    c.W_ro = Real(0.0);
    Real hfus1 = Real(0.0);
    if (c.Aice > p.min_ai) {
        const Real cff = Real(1.0) - brnfr;
        hfus1 = hfus * cff + c.Isst * Cp_w - (cff * Cp_i + brnfr * Cp_w) * c.Tice;
        Qai = c.Qnet_ai;
        const Real Qi2 = b2d * (c.Tice - c.Isst);

        Real Qsur;
        if ((c.Hsno <= eps) && (c.Hmel <= eps)) {
            Qsur = -(Qai - Qi2) / (hfus1 * RhoSW);
        } else if ((c.Hsno <= eps) && (c.Hmel > eps)) {
            Qsur = -(Qai - Qi2) / (hfus1 * Real(1003.1));
        } else {
            Qsur = -(Qai - Qi2) / (hfus * p.SnowWetRho);
        }

        if ((c.Hsno > eps) && (c.Isst >= Real(0.0))) {
            c.Hsno -= c.Aice * std::max(Qsur, Real(0.0)) * dtice;
            c.Hmel += c.Aice * std::max(Qsur, Real(0.0)) * p.SnowWetRho / RhoSW * dtice;
        } else if ((c.Hmel > eps) && (c.Isst <= tfrz)) {
            c.W_ai = std::min(Qsur, Real(0.0));
            c.Hmel += c.Aice * std::min(Qsur, Real(0.0)) * dtice;
        } else if ((c.Hsno <= eps) && (c.Hmel >= eps) && (c.Isst > tfrz)) {
            c.W_ai = std::max(Qsur, Real(0.0));
            c.Hmel += c.Aice * std::max(Qsur, Real(0.0)) * dtice;
        } else if ((c.Hsno < eps) && (c.Hmel < eps) && (c.Isst > tfrz)) {
            c.W_ai = std::max(Qsur, Real(0.0));
            c.Hmel += c.Aice * std::max(Qsur, Real(0.0)) * dtice;
        }

        // Precipitation. ROMS encodes frozen precipitation as rain <= 0.
        if (c.rain <= Real(0.0)) {
            c.Hsno += c.Aice * ws * dtice;
        } else if ((c.Hsno > Real(0.0)) && (c.Hmel == Real(0.0))) {
            c.Hsno = std::max(Real(0.0), c.Hsno - c.Aice * c.rain / p.SnowDryRho);
            c.W_ai -= Real(2.0) * c.Aice * c.rain / p.IceRho;
        } else if ((c.Hsno > Real(0.0)) && (c.Hmel > Real(0.0))) {
            c.Hsno = std::max(Real(0.0), c.Hsno - Real(0.5) * c.Aice * c.rain / p.SnowDryRho);
            c.W_ai -= Real(0.5) * c.Aice * c.rain / p.IceRho;
            c.Hmel += c.Aice * Real(0.5) * c.rain / RhoSW * dtice;
        } else if (c.Hmel > Real(0.0)) {
            c.Hmel += c.Aice * c.rain / RhoSW * dtice;
        } else {
            c.W_ai -= c.Aice * c.rain / p.IceRho;
        }

        // Limit surface water by the smaller of a max limit and the ice
        // thickness; the excess runs off.
        const Real hcap = std::min(p.max_hmelt, c.Hice);
        if (c.Hmel > hcap) {
            c.W_ro = (c.Hmel - hcap) / dtice;
            c.Hmel = hcap;
        }
    }

    // Molecular sublayer under ice. Yaglom and Kader (1974) roughness scales.
    // IceKappa is ROMS's local kappa = 0.4, not REMORA's vonKar = 0.41.
    Real z0 = std::max(z0ii * ice_thick, Real(0.01));
    z0 = std::min(z0, Real(0.1));
    Real zdz0 = c.dztop / z0;
    zdz0 = std::max(zdz0, Real(3.0));
    const Real rno   = utau * Real(0.09) / nu;
    const Real termt = ykf * std::sqrt(rno) * std::pow(prt, Real(0.666667));
    const Real terms = ykf * std::sqrt(rno) * std::pow(prs, Real(0.666667));
    const Real cht = utau / (tpr * (std::log(zdz0) / IceKappa + termt));
    const Real chs = utau / (tpr * (std::log(zdz0) / IceKappa + terms));

    // Temperature and salinity of the molecular sublayer under ice.
    const Real tfz = frln * salt_top;
    c.W_ao = Real(0.0);
    c.W_io = Real(0.0);
    const Real xwai = std::max(Real(0.0), c.W_ai);
    {
        const Real cff = Real(1.0) - brnfr;
        hfus1 = hfus * cff + c.T0mk * Cp_w - (cff * Cp_i + brnfr * Cp_w) * c.Tice;
    }
    if (((temp_top <= tfz) && (c.Qnet_ao > Real(0.0))) ||
        ((temp_top >= tfz) && (c.Qnet_ao < Real(0.0)) && (c.Aice > Real(0.0)))) {
        c.W_ao = c.Qnet_ao / (hfus1 * RhoSW);
    }

    Real xtot;
    if ((c.Aice <= p.min_ai) || (c.Hice <= p.min_hi)) {
        c.S0mk = salt_top;
        c.T0mk = temp_top;
        c.W_ai = Real(0.0);
        xtot = (Real(1.0) - c.Aice) * c.W_ao;
    } else {
        // MK89 version.
        c.W_io = (Qio / RhoSW + Cp_w * cht * (c.T0mk - temp_top)) / hfus1;
        xtot = c.Aice * c.W_io + (Real(1.0) - c.Aice) * c.W_ao;
        // ROMS's note here: this s0mk does not follow from the MK89
        // derivation, but it works well. Four alternatives are commented out
        // in ice_mk.h; this is the one that is live, so it is the one ported.
        c.S0mk = (chs * salt_top + (xwai - c.W_io) * sice) /
                 (chs + xwai + c.W_ro - c.W_io);
        c.S0mk = std::max(c.S0mk, Real(0.0));
        c.S0mk = std::min(c.S0mk, Real(40.0));
        c.T0mk = frln * c.S0mk;
    }

    // Adjust surface heat flux. ROMS's IceCavity branch needs ICESHELF, which is not ported.
    if (c.Aice <= p.min_ai) {
        c.stflx_temp = c.Qnet_ao;
    } else {
        c.stflx_temp = (Real(1.0) - c.Aice) * c.Qnet_ao
                     + c.Aice * Qio
                     - xtot * hfus1 * RhoSW;
    }
    // Back to the ROMS convention (Celsius m/s).
    c.stflx_temp = -c.stflx_temp * RhoCpr;
    c.stflx_temp *= rmask;

    // Adjust surface freshwater flux. On fresh water sice and salt_top are
    // both zero and this term vanishes identically.
    c.stflx_salt += ((c.Aice * (c.W_io - c.W_ai)
                      + (Real(1.0) - c.Aice) * c.W_ao
                      + c.W_fr)) * (salt_top - sice)
                  - c.Aice * (c.W_ro - xwai) * salt_top;
    c.stflx_salt *= rmask;

    // Ice-ocean mass flux.
    c.IOmf = xtot - c.Aice * xwai - c.Aice * c.W_ro + c.W_fr;
    c.IOmf *= rmask;

    // Update ice properties. mi_old tracks the old ice mass so that new
    // production can reduce the average age.
    const Real mi_old = c.Hice;
    Real phi = Real(3.0);
    if (c.W_ao < Real(0.0)) { phi = Real(0.5); }

    c.Hice += dtice * (c.Aice * (c.W_io - c.W_ai)
                       + (Real(1.0) - c.Aice) * c.W_ao + c.W_fr);

    const Real ai_tmp = c.Aice;
    c.Aice += dtice * (Real(1.0) - c.Aice) * (phi * c.W_ao + c.W_fr);
    c.Aice = std::min(c.Aice, p.max_ai);
    if (c.Aice < ai_tmp) {
        c.Hsno *= c.Aice / std::max(ai_tmp, Real(eps));
    }

    // ROMS ICE_CONVSNOW: if the snow base is below sea level, raise it to sea
    // level by converting snow to ice. hstar is weighted by Aice like Hsno
    // and Hice.
    if (p.convsnow) {
        const Real hstar = c.Hsno - (RhoSW - p.IceRho) * c.Hice / p.SnowDryRho;
        if (hstar > Real(0.0)) {
            const Real cff = hstar / RhoSW;
            c.Hsno -= p.IceRho * cff;
            c.Hice += p.SnowDryRho * cff;
        }
    }

    // Age of the sea ice. Any new production reduces the parcel's age.
    if ((c.Iage <= Real(0.0)) && (c.Hice > p.min_hi)) {          // new
        c.Iage = dtice * sec2day;
    } else if ((c.Iage > Real(0.0)) && (c.Hice > p.min_hi)) {    // older
        const Real delta_mi =
            std::min(std::max(c.Hice - mi_old, Real(0.0)) / c.Hice, Real(1.0));
        c.Iage += dtice * sec2day - c.Iage * delta_mi;
    } else {                                                     // melted
        c.Iage = Real(0.0);
    }

    c.Aice *= rmask;
    c.Hice *= rmask;

    // Limit the values.
    c.Aice = std::min(c.Aice, p.max_ai);
    c.Aice = std::max(c.Aice, Real(0.0));
    c.Hice = std::max(c.Hice, Real(0.0));
    c.Hsno = std::max(c.Hsno, Real(0.0));
    c.Hmel = std::max(c.Hmel, Real(0.0));
    c.Tice = std::max(c.Tice, Real(-70.0));
    if (c.Hice <= Real(0.0)) { c.Aice = Real(0.0); }
    if (c.Aice <= Real(0.0)) { c.Hice = Real(0.0); }
}

} // namespace REMORASeaIce

#endif
