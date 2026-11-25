#ifndef REMORA_ICEBULKFLUX_H_
#define REMORA_ICEBULKFLUX_H_

#include <algorithm>
#include <cmath>

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Algorithm.H>

#include <REMORA_SeaIce.H>

namespace REMORASeaIce {

/**
 * The air-ice surface fluxes, ported from the ICE_BULK_FLUXES block of ROMS
 * bulk_flux.F (the section from "Sensible heat flux over ice" to
 * "Qai(i,j)=-Qsh_i-Qlh_i-Qsw_i+Qlw_i").
 *
 * This is the other half of the coupling the ice model needs. REMORA's
 * bulk_fluxes() already produces the open-water flux -- ROMS's Qnet_ao is
 * exactly -stflux(itemp)/Hscale (bulk_flux.F:1283) -- so what has to be added
 * is the flux over the ICE surface, plus the two coefficients with which
 * ice_mk.h solves for the ice surface temperature implicitly:
 *
 *   Qai       net heat flux from ice to atmosphere [W/m2]
 *   Coef_Ice  d(flux)/d(Tice), linearised            [(W/m2)/K] -> Fi(icQcon)
 *   RHS_Ice   the constant part of the same          [W/m2]     -> Fi(icQrhs)
 *
 * The thermodynamics then forms Isst = (Qrhs/Qcon) - 273.15, which is why
 * everything here is in Kelvin.
 *
 * Three constants in this block are NOT the ones a reader might expect, and
 * the port keeps ROMS's:
 *   - `blk_Cpa` (1004.67, Businger 1982) is used here, not `ice.in`'s
 *     `spec_heat_air` (1004.0). Together with d1/d2i being dead in ice_mk.h
 *     that makes `spec_heat_air` unreachable in this configuration.
 *   - `emmiss` (0.97), the OCEAN emissivity from mod_scalars, is used for the
 *     ice longwave -- not `ice_emiss`. Both are 0.97, so there is no
 *     numerical consequence today, but they are different inputs.
 *   - Ch and Ce are a hard-wired 1.75e-3, attributed in ROMS to "Tale Bakken
 *     Ulfsby, 2025"; the COARE-iteration forms are commented out and dated
 *     02/2024. They are NOT derived from the stress iteration.
 */
constexpr amrex::Real IceTransferCoeff = amrex::Real(1.75e-3);

struct IceAirFluxIn {
    amrex::Real rhoAir;    //!< air density [kg/m3]
    amrex::Real delW;      //!< effective wind speed incl. gustiness [m/s]
    amrex::Real Q;         //!< air specific humidity [kg/kg]
    amrex::Real Qsea;      //!< saturation specific humidity at the water [kg/kg]
    amrex::Real TairC;     //!< air temperature [C]
    amrex::Real Hlv;       //!< latent heat of vaporisation [J/kg]
    amrex::Real Pair;      //!< air pressure [mb]
    amrex::Real rain;      //!< precipitation [kg/m2/s]; ROMS reads rain < 0 as frozen
    int  lw_mode;          //!< 0 computed (ROMS LONGWAVE), 1 downward given (LONGWAVE_OUT), 2 net given
    amrex::Real lw_in;     //!< the longwave given as input [W/m2]: downward in mode 1, net in mode 2
    amrex::Real cloud;     //!< cloud fraction, read in mode 0
    amrex::Real srflx_ice; //!< net shortwave reaching the ice surface [W/m2], before the ice attenuates it;
                           //!< not REMORA's vec_srflx_ice, which is the shortwave reaching the water
    amrex::Real albedo_ice;//!< ice/snow albedo
    amrex::Real wet_bulb;  //!< wet-bulb factor, from REMORA's bulk_fluxes
    amrex::Real temp_top;  //!< surface water temperature [C]
    amrex::Real Aice;      //!< ice concentration [fraction]
    amrex::Real Hmel;      //!< melt-water thickness [m]
    amrex::Real Isst;      //!< ice surface temperature from the previous step [C]
    amrex::Real blk_ZT;    //!< air-temperature measurement height [m]
};

struct IceAirFluxOut {
    amrex::Real Qai;       //!< net heat flux, ice to atmosphere [W/m2]
    amrex::Real Coef_Ice;  //!< [(W/m2)/K], seeds Fi(icQcon)
    amrex::Real RHS_Ice;   //!< [W/m2], seeds Fi(icQrhs)
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
ice_air_flux (const IceAirFluxIn& in,
              const SeaIceParameters& p,
              IceAirFluxOut& out)
{
    using amrex::Real;

    // ROMS's own values for this block, not ice.in's; see the note above.
    constexpr Real blk_Cpa = Real(1004.67);
    constexpr Real emmiss  = Real(0.97);
    // ROMS mod_scalars.F. REMORA's own StefBo is the CODATA value, 6.6e-5 larger.
    constexpr Real StefBo  = Real(5.67e-8);
    const Real Ch = IceTransferCoeff;
    const Real Ce = IceTransferCoeff;

    const Real TairK = in.TairC + Real(273.16);
    const Real TiceK = (in.Aice > p.min_ai) ? (in.Isst     + Real(273.16))
                                            : (in.temp_top + Real(273.16));

    Real Coef = Real(0.0);
    Real RHS  = Real(0.0);

    // Sensible heat flux over ice.
    Real cff = in.rhoAir * blk_Cpa * Ch * in.delW;              // (W/m2)/K
    Real Qsh_i = cff * (TairK + Real(0.0098) * in.blk_ZT - TiceK);
    Coef += cff;
    RHS  += cff * (TairK + Real(0.0098) * in.blk_ZT);

    // Sensible heat carried by warm rain or cold snow reaching the ice, using
    // the same wet-bulb factor as the open-water sensible flux.
    cff = in.wet_bulb * Real(2093.0);
    Qsh_i += std::abs(in.rain) * cff
           * (TiceK - in.TairC + (in.Qsea - in.Q) * in.Hlv / blk_Cpa);
    Coef -= std::abs(in.rain) * cff;
    // TairK here rather than TairC, because Isst is solved in Kelvin.
    RHS  += std::abs(in.rain) * cff
          * (-TairK + (in.Qsea - in.Q) * in.Hlv / blk_Cpa);

    // Latent heat flux over ice. Saturation specific humidity over ice from
    // Parkinson and Washington (1979).
    constexpr Real le_i = Real(2.834e6);
    const Real slp = in.Pair * Real(100.0);                     // mb -> Pa
    cff = Real(611.0) * std::pow(Real(10.0),
              Real(9.5) * (TiceK - Real(273.16)) / (TiceK - Real(7.66)));
    const Real Qsati = Real(0.62197) * cff
                     / (slp - (Real(1.0) - Real(0.62197)) * cff);
    const Real dq_i = in.Q - Qsati;
    Real Qlh_i = in.rhoAir * Ce * in.delW * le_i * dq_i;
    // Rain freezing on the surface releases latent heat, but only where there
    // is no melt water: with a pond present ROMS assumes the rain just pools.
    if ((in.rain > Real(0.0)) && (in.Hmel == Real(0.0))) {
        Qlh_i += Thermo::hfus * in.rain;
    }
    RHS += Qlh_i;

    // Shortwave to the ice. ROMS's note: it is an open question what fraction
    // of this should drive the surface temperature rather than the interior;
    // all of it does, here as there.
    const Real Qsw_i = (Real(1.0) - in.albedo_ice) * in.srflx_ice;
    RHS += Qsw_i;

    // Net longwave from the ice, by the same route the ocean's longwave takes.
    // Qlw_i is the flux from ice to atmosphere.
    Real Qlw_i;
    if (in.lw_mode == 0) {
        // Computed from the air temperature, humidity and cloud (ROMS LONGWAVE)
        cff = in.Q / (Real(1.0) + in.Q);
        const Real vap_p_i = slp * cff / (Real(0.62197) + cff);
        const Real cff1 = Real(0.39) - Real(0.005) * std::sqrt(vap_p_i);
        const Real cff2 = Real(1.0) - Real(0.65) * in.cloud;
        const Real cff3 = StefBo * emmiss * (TairK * TairK * TairK);
        const Real cff4 = cff3 * TairK;
        Qlw_i = cff4 * cff1 * cff2 + Real(4.0) * cff3 * (TiceK - TairK);
        Coef += Real(4.0) * cff3;
        RHS  -= cff4 * (cff1 * cff2 - Real(4.0));
    } else if (in.lw_mode == 1) {
        // Downward longwave given (ROMS LONGWAVE_OUT)
        const Real cff3 = StefBo * emmiss * (TiceK * TiceK * TiceK);
        const Real cff4 = cff3 * TiceK;
        Qlw_i = -(in.lw_in - cff4);
        Coef += Real(4.0) * cff3;
        RHS  += in.lw_in + Real(3.0) * cff4;
    } else {
        // Net longwave given (ROMS with neither option)
        const Real cff3 = StefBo * emmiss * (TairK * TairK * TairK);
        Qlw_i = in.lw_in;
        Coef += Real(4.0) * cff3;
        RHS  -= in.lw_in + Real(4.0) * cff3 * (TiceK - Real(2.0) * TairK);
    }

    out.Qai      = -Qsh_i - Qlh_i - Qsw_i + Qlw_i;
    out.Coef_Ice = Coef;
    out.RHS_Ice  = RHS;
}

/** Open-water albedo, the constant alb_w of ROMS ice_albedo.F */
constexpr amrex::Real alb_w = amrex::Real(0.06);

/**
 * Ice and snow albedo, ported from the DEFAULT formulation of ROMS
 * ice_albedo.F (neither ICE_ALB_EC92 nor ALBEDO_CSIM).
 *
 * **ROMS's melting branch is unreachable, and this reproduces that by
 * default.** ice_albedo.F takes `sfc_temp = Fi(i,j,icIsst)` and tests
 * `(sfc_temp - 273.16) > -1.0` for the wet albedo, but Isst is in CELSIUS:
 * ice_mk.h writes `(Qrhs/Qcon) - 273.15` and clamps it to [-45, 0], and
 * the only read that converts is bulk_flux.F's `TiceK = Fi(icIsst) + 273.16`.
 * So `sfc_temp - 273.16` lies in [-318, -273] and is always below -1: the
 * albedo is always the dry value, 0.85 over snow and 0.65 over bare ice.
 *
 * `melt_fix = true` applies the evidently intended Celsius comparison instead.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
ice_albedo (const amrex::Real Aice,
            const amrex::Real Hsno,
            const amrex::Real Isst,
            const SeaIceParameters& p,
            const bool melt_fix)
{
    using amrex::Real;
    constexpr Real alb_i_dry = Real(0.65);
    constexpr Real alb_i_wet = Real(0.60);
    constexpr Real alb_s_dry = Real(0.85);
    constexpr Real alb_s_wet = Real(0.72);

    if (Aice <= p.min_ai) { return alb_w; }

    const Real cff1 = alb_s_wet - alb_s_dry;
    const Real cff2 = alb_i_wet - alb_i_dry;
    // ROMS: sfc_temp - 273.16 > -1, with sfc_temp in Celsius. melt_fix drops
    // the spurious conversion and compares the Celsius value directly.
    const bool melting = melt_fix ? (Isst > Real(-1.0))
                                  : ((Isst - Real(273.16)) > Real(-1.0));
    // The reference temperature in the slope term is 272.16 K either way; with
    // the fix it becomes -1 C, which is the same physical point.
    const Real tref = melt_fix ? Real(-1.0) : Real(272.16);

    if (Hsno > Real(0.0)) {
        return melting ? (cff1 * (Isst - tref) + alb_s_dry) : alb_s_dry;
    }
    return melting ? (cff2 * (Isst - tref) + alb_i_dry) : alb_i_dry;
}

/**
 * Ice and snow albedo of Ebert and Curry (1992), ROMS ice_albedo.F under
 * ICE_ALB_EC92. The thicknesses are per unit ice area, formed as ROMS forms
 * them here, with 0.001 added to the concentration.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
ice_albedo_ec92 (const amrex::Real Aice,
                 const amrex::Real Hice,
                 const amrex::Real Hsno,
                 const amrex::Real Hmel,
                 const SeaIceParameters& p)
{
    using amrex::Real;
    if (Aice <= p.min_ai) { return alb_w; }

    const Real cff = Real(1.0) / (Aice + Real(0.001));
    // The fit goes negative for very thin ice, so ROMS floors the thickness.
    const Real hi = amrex::max(cff * Hice, Real(0.01));
    const Real hs = cff * Hsno;

    Real alb_ice;
    if (hi >= Real(2.0)) {
        alb_ice = Real(0.561632);
    } else if (hi >= Real(1.0)) {
        alb_ice = Real(0.07616) * hi + Real(0.414492);
    } else {
        alb_ice = Real(0.082409) * std::log(hi) + Real(0.485472);
    }

    Real alb_snow;
    if (Hmel > Real(0.0)) {
        alb_snow = (hs >= Real(0.1)) ? Real(0.701009)
                                     : alb_ice + (Real(0.701009) - alb_ice) * hs / Real(0.1);
    } else {
        alb_snow = Real(0.83);
    }

    if (Hsno > Real(0.0)) { return alb_snow; }
    if (Hmel > Real(0.02)) { return Real(0.42); }
    return alb_ice;
}

} // namespace REMORASeaIce

#endif
