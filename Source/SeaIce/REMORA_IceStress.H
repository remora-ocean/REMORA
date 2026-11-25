#ifndef REMORA_ICESTRESS_H_
#define REMORA_ICESTRESS_H_

#include <algorithm>
#include <cmath>

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Algorithm.H>

#include <REMORA_SeaIce.H>

namespace REMORASeaIce {

/**
 * Per-point pieces of ROMS ice_vbc.F.
 *
 * ice_vbc's only stencil is a 2-point average -- (i-1,i) at u points and
 * (j-1,j) at v points -- so these take the already-averaged quantities and
 * the caller does the averaging.
 *
 * Two constants differ from ice_mk.h's despite sharing a name, and the port
 * keeps ROMS's:
 *   - `kappa` is 0.4 here as there (NOT REMORA's vonKar 0.41).
 *   - `z0ii` is **0.01 m** here and **0.02 m** in ice_mk.h. Same name, two
 *     values, in two files of the same module. Not a typo to fix: each is
 *     used only in its own routine.
 *
 * Deliberately absent: `rhoO`, which ice_vbc.F computes in all three of its
 * loops and never reads.
 */
namespace VBC {
    constexpr amrex::Real kappa = amrex::Real(0.4);     //!< ROMS ice_vbc.F, local
    constexpr amrex::Real z0ii  = amrex::Real(0.01);    //!< ice_vbc.F -- NOT ice_mk.h's 0.02
    constexpr amrex::Real eps   = amrex::Real(1.0e-20);
}

/**
 * The rho-point coefficients: friction velocity, ice-ocean transfer
 * coefficient, and the air-ice form-drag weighting.
 *
 * `IOmt_prev` is the PREVIOUS step's chuiw, which is what ROMS reads out of
 * Fi(icIOmt) here before overwriting it at the end of the routine.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
ice_vbc_coef (const amrex::Real IOvs,
              const amrex::Real Hice,
              const amrex::Real Aice,
              const amrex::Real dztop,
              const amrex::Real IOmt_prev,
              const SeaIceParameters& p,
              amrex::Real& utau,
              amrex::Real& chuiw,
              amrex::Real& chuai)
{
    using amrex::Real;
    Real spd = amrex::max(IOvs, Real(0.15));
    if (p.mom_bulk) {
        // ROMS ICE_MOM_BULK. The air-ice stress carries no form-drag factor
        // in this branch; a factor of exactly one leaves ice_vbc_stress's
        // product unchanged.
        utau  = spd;
        chuiw = p.Cd_io * spd;
        chuai = Real(1.0);
        return;
    }
    const Real thic = Hice / amrex::max(Aice, amrex::max(p.min_ai, Real(VBC::eps)));
    Real z0 = amrex::max(VBC::z0ii * thic, Real(0.01));
    z0 = amrex::min(z0, Real(0.1));
    Real zdz0 = dztop / z0;
    if (zdz0 < Real(6.0)) { zdz0 = Real(6.0); }
    utau = std::sqrt(IOmt_prev * spd);
    utau = amrex::max(utau, Real(1.0e-4));
    chuiw = VBC::kappa * utau / std::log(zdz0);
    chuai = Real(1.0) - std::cos(Real(1.0) * Real(3.14159265358979323846)
                                 * amrex::min(thic + Real(0.05), Real(1.0)));
}

/**
 * The stress on the WATER at one velocity point, and the air-ice stress on
 * the ice. ROMS:
 *
 *   sustr = aix*chux*(Uice - Uavg) + (1 - aix)*sustr_ao
 *   AIus  = 0.5*aix*chuax*(sustr_ai(i) + sustr_ai(i-1))/IceRho
 *
 * At full cover the wind stress on the water is REPLACED by the ice-water
 * drag. ROMS's IceCavity branch is unconditionally false without ICESHELF, so
 * only this assignment exists.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
ice_vbc_stress (const amrex::Real ai_avg,    //!< concentration averaged to the face
                const amrex::Real chu_avg,   //!< chuiw averaged to the face
                const amrex::Real chua_avg,  //!< chuai averaged to the face
                const amrex::Real vice,      //!< ice velocity at the face
                const amrex::Real vavg,      //!< mixed-layer velocity at the face
                const amrex::Real str_ao,    //!< open-water wind stress at the face
                const amrex::Real str_ai_hi, //!< air-ice stress at (i,j)
                const amrex::Real str_ai_lo, //!< air-ice stress at (i-1,j) or (i,j-1)
                const SeaIceParameters& p,
                amrex::Real& str_water,
                amrex::Real& str_on_ice)
{
    using amrex::Real;
    str_on_ice = Real(0.5) * ai_avg * chua_avg
               * (str_ai_hi + str_ai_lo) / p.IceRho;
    str_water  = ai_avg * chu_avg * (vice - vavg)
               + (Real(1.0) - ai_avg) * str_ao;
}

} // namespace REMORASeaIce

#endif
