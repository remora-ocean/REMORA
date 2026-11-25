#ifndef REMORA_ICEEVP_H_
#define REMORA_ICEEVP_H_

#include <algorithm>
#include <cmath>

#include <AMReX_REAL.H>
#include <AMReX_Algorithm.H>
#include <AMReX_GpuQualifiers.H>

#include <REMORA_SeaIce.H>

namespace REMORASeaIce {

/**
 * Elastic-viscous-plastic rheology, ported from ROMS ice_evp.F and
 * ice_evp_sig.F (Hunke and Dukowicz 1997, as Budgell 2005 implements it).
 *
 * ROMS computes the strain rates in both routines, from the same expressions; here they are
 * in ice_strain_rates(). The stencil does not depend on the land mask, unlike ice_elastic.F's.
 *
 * Components of the internal ice stress:
 * ROMS keeps them in `Si` and indexes their two time levels separately from the rest of it.
 * Here the stress has its own field and two MultiFabs.
 */
namespace IceSig {
    enum { xx = 0, yy, xy, NumComps };
}

namespace EVP {
    constexpr amrex::Real epso = amrex::Real(1.0e-12);  //!< delta floor, ice_evp.F
    constexpr amrex::Real E0   = amrex::Real(0.25);     //!< ice_evp_sig.F
}

/**
 * The three strain-rate components at a rho point.
 *
 * ROMS, with `lieol` the current EVP sub-level:
 *
 *   eps_xx = (Uevp(i+1,j) - Uevp(i,j)) * pm
 *   eps_yy = (Vevp(i,j+1) - Vevp(i,j)) * pn
 *   epx    = 0.25*(Vevp(i+1,j+1)+Vevp(i+1,j)-Vevp(i-1,j+1)-Vevp(i-1,j)) * pm
 *   epy    = 0.25*(Uevp(i+1,j+1)+Uevp(i,j+1)-Uevp(i+1,j-1)-Uevp(i,j-1)) * pn
 *   eps_xy = 0.5*(epx + epy)
 *
 * `epx` is the centred derivative dv/dx at the rho point: the factor 0.25 is
 * one half for averaging the two v points onto the rho row and one half for
 * the two-cell span from i-1 to i+1. `epy` is du/dy in the same way, so
 * eps_xy is the usual 0.5*(du/dy + dv/dx).
 */
struct IceStrain { amrex::Real xx, yy, xy; };

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
IceStrain
ice_strain_rates (const amrex::Real u_i_j,    const amrex::Real u_ip1_j,
                  const amrex::Real u_i_jp1,  const amrex::Real u_ip1_jp1,
                  const amrex::Real u_i_jm1,  const amrex::Real u_ip1_jm1,
                  const amrex::Real v_i_j,    const amrex::Real v_i_jp1,
                  const amrex::Real v_ip1_j,  const amrex::Real v_ip1_jp1,
                  const amrex::Real v_im1_j,  const amrex::Real v_im1_jp1,
                  const amrex::Real pm,       const amrex::Real pn)
{
    using amrex::Real;
    IceStrain e;
    e.xx = (u_ip1_j - u_i_j) * pm;
    e.yy = (v_i_jp1 - v_i_j) * pn;
    const Real epx = Real(0.25) * (v_ip1_jp1 + v_ip1_j - v_im1_jp1 - v_im1_j) * pm;
    const Real epy = Real(0.25) * (u_ip1_jp1 + u_i_jp1 - u_ip1_jm1 - u_i_jm1) * pn;
    e.xy = Real(0.5) * (epx + epy);
    return e;
}

/**
 * Ice strength and the bulk and shear viscosities (ROMS ice_evp.F).
 *
 * **ROMS computes a ceiling it does not use.** `zmax = 2.5e8*Pice` is
 * evaluated and never read; the clamp is `zetaMax`, whose ice.in default is
 * 1e14. Wherever `delta` sits at its 1e-12 floor -- near-rigid ice -- that
 * leaves the viscosity 20x to 2000x stiffer than `zmax` would allow. The
 * default here reproduces ROMS; `zmax_clamp = true` applies the ceiling ROMS
 * wrote down.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
ice_evp_viscosity (const IceStrain& e,
                   const amrex::Real Hice,
                   const amrex::Real Aice,
                   const SeaIceParameters& p,
                   const bool zmax_clamp,
                   amrex::Real& Pice,
                   amrex::Real& Bvis,
                   amrex::Real& Svis,
                   const amrex::Real Pgrd = amrex::Real(0.0))
{
    using amrex::Real;
    const Real e2r = Real(1.0) / p.ellip_sq;

    const Real eone  = e.xx + e.yy;
    const Real etwos = (e.xx - e.yy) * (e.xx - e.yy)
                     + Real(4.0) * e.xy * e.xy;
    Real delta = std::abs(eone * eone + e2r * etwos);
    delta = amrex::max(std::sqrt(delta), Real(EVP::epso));

    // ROMS ICE_STRENGTH_QUAD takes the strength from a gridded parameter and
    // the square of the thickness. ROMS never sets that parameter, so with
    // the option on the strength is zero and the viscosity sits at zetaMin.
    Pice = (p.strength_quad ? Pgrd * Hice * Hice : p.Pstar * Hice)
         * std::exp(-p.Astrength * (Real(1.0) - Aice));
    Bvis = Pice / (Real(2.0) * delta);
    const Real zmax = Real(2.5e8) * Pice;      // ROMS computes this and never reads it
    if (zmax_clamp) { Bvis = amrex::min(Bvis, zmax); }
    Bvis = amrex::min(Bvis, p.zetaMax);
    Bvis = amrex::max(Bvis, p.zetaMin);
    Svis = e2r * Bvis;
}

/**
 * Time-step the internal ice stress (ROMS ice_evp_sig.F).
 *
 * Above a thickness of 0.01 m the elastic update runs; at or below it ROMS
 * falls back to the pure viscous-plastic stress. `dtevp` is `dtice/nEVP`.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
ice_evp_stress (const IceStrain& e,
                const amrex::Real Hice,
                const amrex::Real Pice,
                const amrex::Real Bvis,
                const amrex::Real Svis,
                const amrex::Real sxx_old,
                const amrex::Real syy_old,
                const amrex::Real sxy_old,
                const amrex::Real pm,
                const amrex::Real dtevp,
                const SeaIceParameters& p,
                amrex::Real& sxx_new,
                amrex::Real& syy_new,
                amrex::Real& sxy_new)
{
    using amrex::Real;
    if (Hice > Real(0.01)) {
        const Real pmdt = pm * dtevp;
        const Real E   = Real(2.0) * EVP::E0 * p.IceRho * Hice / (pmdt * pmdt);
        const Real ep  = E * Pice / (Real(4.0) * Bvis + Real(1.0e-8));
        const Real ee  = E / (Real(2.0) * Svis + Real(1.0e-8));
        const Real ees = E * (Svis - Bvis)
                       / (Real(4.0) * Svis * Bvis + Real(1.0e-8));
        const Real cff   = Real(1.0) / dtevp;
        const Real alfa  = cff + ee + ees;
        const Real beta  = ees;
        const Real gamma = Real(1.0) / dtevp + ee;
        const Real f1 = E * e.xx + cff * sxx_old - ep;
        const Real f2 = E * e.yy + cff * syy_old - ep;
        const Real f3 = E * e.xy + cff * sxy_old;
        const Real f4 = Real(1.0) / (alfa * alfa - beta * beta);
        sxx_new = f4 * (alfa * f1 - beta * f2);
        sxy_new = f3 / gamma;
        syy_new = f4 * (alfa * f2 - beta * f1);
    } else {
        const Real div = e.xx + e.yy;
        sxx_new = Real(2.0) * Svis * e.xx + (Bvis - Svis) * div
                - Real(0.5) * Pice;
        syy_new = Real(2.0) * Svis * e.yy + (Bvis - Svis) * div
                - Real(0.5) * Pice;
        sxy_new = Real(2.0) * Svis * e.xy;
    }
}

} // namespace REMORASeaIce

#endif
