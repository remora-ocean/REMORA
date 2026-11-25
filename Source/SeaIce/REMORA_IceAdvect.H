#ifndef REMORA_ICEADVECT_H_
#define REMORA_ICEADVECT_H_

#include <algorithm>
#include <cmath>

#include <AMReX_REAL.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>

#include <REMORA_SeaIce.H>

namespace REMORASeaIce {

/**
 * MPDATA advection of the ice tracers, ported from ROMS ice_smolar.h
 * (Smolarkiewicz's multidimensional positive-definite scheme).
 *
 * Two passes: a donor-cell upwind step, then an antidiffusive correction
 * computed from a pseudo-velocity. ROMS runs it once per advected field (Aice, Hice, Hsno,
 * Hmel, the enthalpy and Hage, from which Iage is recovered), and `ice_limit.F` does the
 * clipping afterwards.
 *
 * **The scheme is positive-definite, not sign-agnostic.** The `rateu` and
 * `rateyiu` denominators are `MAX(epsil, <sum of the field>)`, which is only
 * a guard if the field is non-negative. That is why `ice_advect.F` negates
 * the enthalpy before advecting it and negates it back afterwards -- the
 * enthalpy is the one advected field that is negative.
 *
 * ROMS's `ICE_UPWIND` leaves out the antidiffusive pass (`remora.ice.advect_upwind`).
 */
namespace MPDATA {
    constexpr amrex::Real epsil = amrex::Real(1.0e-15);   //!< ice_smolar.h
}

/** Donor-cell upwind flux across a u face. `on_u` is the face length. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
mpdata_upwind_u (const amrex::Real ui,
                 const amrex::Real f_im1, const amrex::Real f_i,
                 const amrex::Real on_u)
{
    using amrex::Real;
    const Real cff1 = amrex::max(Real(0.0), ui);
    const Real cff2 = amrex::min(Real(0.0), ui);
    return on_u * (cff1 * f_im1 + cff2 * f_i);
}

/** Donor-cell upwind flux across a v face. `om_v` is the face length. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
mpdata_upwind_v (const amrex::Real vi,
                 const amrex::Real f_jm1, const amrex::Real f_j,
                 const amrex::Real om_v)
{
    using amrex::Real;
    const Real cff1 = amrex::max(Real(0.0), vi);
    const Real cff2 = amrex::min(Real(0.0), vi);
    return om_v * (cff1 * f_jm1 + cff2 * f_j);
}

/**
 * The antidiffusive flux across a u face.
 *
 * `aif` is the field after the upwind pass, `FE` the masked j-differences of
 * it. ROMS:
 *
 *   rateu   = (aif(i) - aif(i-1)) / MAX(epsil, aif(i) + aif(i-1))
 *   rateyiu = sum of four FE / MAX(epsil, <the same cells plus their FE
 *                                   divergence>)
 *   Cu      = 0.5*dt*(pm(i)+pm(i-1))*ui
 *   Cu_crss = 0.5*dt*0.0625*(four pn)*(four vi)
 *   uspeed  = rateu*(|ui| - Cu*ui) - rateyiu*Cu_crss*ui
 *
 * and then upwinds `aif` with `uspeed` in place of the velocity.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
mpdata_anti_u (const int i, const int j,
               const amrex::Array4<amrex::Real const>& aif,
               const amrex::Array4<amrex::Real const>& FE,
               const amrex::Array4<amrex::Real const>& ui,
               const amrex::Array4<amrex::Real const>& vi,
               const amrex::Array4<amrex::Real const>& pm,
               const amrex::Array4<amrex::Real const>& pn,
               const amrex::Real on_u,
               const amrex::Real dtice)
{
    using amrex::Real;
    const Real a_i = aif(i,j,0), a_m = aif(i-1,j,0);
    const Real rateu = (a_i - a_m) / amrex::max(Real(MPDATA::epsil), a_i + a_m);
    const Real rateyiu =
        (FE(i,j+1,0) + FE(i,j,0) + FE(i-1,j+1,0) + FE(i-1,j,0))
        / amrex::max(Real(MPDATA::epsil),
                     a_i + FE(i  ,j+1,0) - FE(i  ,j,0)
                   + a_m + FE(i-1,j+1,0) - FE(i-1,j,0));
    const Real Cu = Real(0.5) * dtice * (pm(i,j,0) + pm(i-1,j,0)) * ui(i,j,0);
    const Real Cu_crss = Real(0.5) * dtice * Real(0.0625)
        * (pn(i-1,j+1,0) + pn(i,j+1,0) + pn(i-1,j-1,0) + pn(i,j-1,0))
        * (vi(i-1,j+1,0) + vi(i,j+1,0) + vi(i-1,j  ,0) + vi(i,j  ,0));
    const Real uspeed = rateu * (std::abs(ui(i,j,0)) - Cu * ui(i,j,0))
                      - rateyiu * Cu_crss * ui(i,j,0);
    const Real cff1 = amrex::max(Real(0.0), uspeed);
    const Real cff2 = amrex::min(Real(0.0), uspeed);
    return on_u * (cff1 * a_m + cff2 * a_i);
}

/** The antidiffusive flux across a v face; the mirror of mpdata_anti_u. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
mpdata_anti_v (const int i, const int j,
               const amrex::Array4<amrex::Real const>& aif,
               const amrex::Array4<amrex::Real const>& FX,
               const amrex::Array4<amrex::Real const>& ui,
               const amrex::Array4<amrex::Real const>& vi,
               const amrex::Array4<amrex::Real const>& pm,
               const amrex::Array4<amrex::Real const>& pn,
               const amrex::Real om_v,
               const amrex::Real dtice)
{
    using amrex::Real;
    const Real a_j = aif(i,j,0), a_m = aif(i,j-1,0);
    const Real ratev = (a_j - a_m) / amrex::max(Real(MPDATA::epsil), a_j + a_m);
    const Real ratexiv =
        (FX(i+1,j,0) + FX(i,j,0) + FX(i+1,j-1,0) + FX(i,j-1,0))
        / amrex::max(Real(MPDATA::epsil),
                     a_j + FX(i+1,j  ,0) - FX(i,j  ,0)
                   + a_m + FX(i+1,j-1,0) - FX(i,j-1,0));
    const Real Cu = Real(0.5) * dtice * (pn(i,j,0) + pn(i,j-1,0)) * vi(i,j,0);
    const Real Cu_crss = Real(0.5) * dtice * Real(0.0625)
        * (pm(i+1,j,0) + pm(i+1,j-1,0) + pm(i-1,j,0) + pm(i-1,j-1,0))
        * (ui(i,j  ,0) + ui(i+1,j  ,0) + ui(i,j-1,0) + ui(i+1,j-1,0));
    const Real vspeed = ratev * (std::abs(vi(i,j,0)) - Cu * vi(i,j,0))
                      - ratexiv * Cu_crss * vi(i,j,0);
    const Real cff1 = amrex::max(Real(0.0), vspeed);
    const Real cff2 = amrex::min(Real(0.0), vspeed);
    return om_v * (cff1 * a_m + cff2 * a_j);
}

/**
 * The state limiter, ported from ROMS ice_limit.F.
 *
 * **This is NOT the same limiter as the one inline at the end of ice_mk.h**,
 * and ROMS applies both -- the inline one inside the thermodynamics, then
 * this one from seaice.F, after the advection and again after the
 * thermodynamics. Three differences, all of which matter:
 *
 *   - this one also clamps Tice from ABOVE, at 0 C; ice_mk.h only floors it
 *     at -70
 *   - this one floors Iage at 0; ice_mk.h does not touch Iage
 *   - when the ice vanishes, this one zeroes FIVE fields -- Aice, Hice, Hmel,
 *     Hsno and Iage -- and triggers on EITHER Aice <= 0 OR Hice <= 0, where
 *     ice_mk.h zeroes only the other of the two
 *
 * So the melt water, the snow and the age survive a vanishing of the ice
 * through the thermodynamics and are cleared by this. Writing one limiter for
 * both would be wrong in both places.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
ice_limit_state (amrex::Real& Aice, amrex::Real& Hice,
                 amrex::Real& Hsno, amrex::Real& Hmel,
                 amrex::Real& Tice, amrex::Real& Iage,
                 const SeaIceParameters& p)
{
    using amrex::Real;
    Aice = amrex::min(Aice, p.max_ai);
    Aice = amrex::max(Aice, Real(0.0));
    Hice = amrex::max(Hice, Real(0.0));
    Hsno = amrex::max(Hsno, Real(0.0));
    Hmel = amrex::max(Hmel, Real(0.0));
    Tice = amrex::max(Tice, Real(-70.0));
    Tice = amrex::min(Tice, Real(0.0));
    Iage = amrex::max(Iage, Real(0.0));
    if ((Aice <= Real(0.0)) || (Hice <= Real(0.0))) {
        Aice = Real(0.0);
        Hice = Real(0.0);
        Hmel = Real(0.0);
        Hsno = Real(0.0);
        Iage = Real(0.0);
    }
}

} // namespace REMORASeaIce

#endif
