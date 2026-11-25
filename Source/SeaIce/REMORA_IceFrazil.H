#ifndef REMORA_ICEFRAZIL_H_
#define REMORA_ICEFRAZIL_H_

#include <algorithm>
#include <cmath>

#include <AMReX_REAL.H>
#include <AMReX_Algorithm.H>
#include <AMReX_GpuQualifiers.H>

#include <REMORA_SeaIce.H>

namespace REMORASeaIce {

/**
 * Frazil ice, ported from ROMS ice_frazil.F.
 *
 * The one part of the ice model that reaches into the 3D ocean: it walks the
 * water column, forms ice wherever the water is supercooled and warms that
 * water back toward the freezing point, then melts the accumulated frazil
 * again in any layer that is above freezing. The accumulated rate goes to
 * Fi(icW_fr), which the thermodynamics reads.
 *
 * ROMS never calls ice_frazil.F: its only call is inside the file itself. Without frazil a
 * freshwater column, which freezes near 0 C at every depth, supercools without bound. It is
 * off by default, which is what ROMS runs; `remora.ice.do_frazil = true` turns it on. Having
 * no ROMS ordering to match, it runs after the ocean step, and the rate it leaves is read by
 * the next step's thermodynamics.
 *
 * Two things in ROMS worth knowing, both reproduced as written:
 *
 * 1. **The freezing point has no pressure dependence.** ROMS's
 *    `freezing_point(S, Z)` takes a depth and ignores it: the live line is
 *    `FP = -0.0543*S`, with both the Gill (1982) and Steele et al. (1989)
 *    forms -- which carry `+0.000753*Z` and `+0.000759*Z` -- commented out.
 *    At 100 m the omitted term is -0.076 C, which is not small for a
 *    supercooling process. `pressure_freezing = true` restores the Steele
 *    form; the default reproduces ROMS.
 *
 * 2. **The two branches use different reference densities.** The freezing
 *    branch weights by `(1000.0 + rho)` and the melting branch by
 *    `(rho0 + rho)`. They coincide at rho0 = 1000; at rho0 = 1025 freezing and
 *    melting are weighted 2.5 % apart.
 */
namespace Frazil {
    constexpr amrex::Real Lhat = amrex::Real(79.2);   //!< ice_frazil.F
    constexpr amrex::Real r    = amrex::Real(0.5);    //!< ice_frazil.F
    constexpr amrex::Real pfac = amrex::Real(1.0);    //!< ROMS's "adjusting factor"
}

/**
 * Freezing point of sea water.
 *
 * `z` is negative downward and is IGNORED unless `pressure` is set, which is
 * ROMS's behaviour -- see note 1 above.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
freezing_point (const amrex::Real S, const amrex::Real z, const bool pressure)
{
    using amrex::Real;
    // Steele et al. (1989), the form ROMS has commented out.
    return pressure ? (Real(-0.0543) * S + Real(0.000759) * z)
                    : (Real(-0.0543) * S);
}

/**
 * One level of the frazil calculation.
 *
 * `w_fr` accumulates across the column (the caller walks k from the bottom
 * up, as ROMS does) and `temp` is updated in place. `at_surface` is true at
 * k = N, where ROMS suppresses a negative gamma.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
ice_frazil_level (amrex::Real& temp, const amrex::Real salt,
                  const amrex::Real z_r, const amrex::Real Hz,
                  const amrex::Real rho, const amrex::Real rho0,
                  const bool at_surface,
                  const SeaIceParameters& p,
                  const bool pressure_freezing,
                  amrex::Real& w_fr)
{
    using amrex::Real;
    const Real orhoi = Real(1.0) / p.IceRho;
    const Real t_fr = freezing_point(salt, z_r, pressure_freezing);

    if (temp < t_fr) {
        Real gamma_k = Frazil::pfac * (t_fr - temp)
            / (Frazil::Lhat + temp * (Real(1.0) - Frazil::r)
               + Real(0.0543) * salt);
        if ((gamma_k < Real(0.0)) && at_surface) { gamma_k = Real(0.0); }
        // NOTE ROMS writes 1000.0 here and rho0 in the melting branch below.
        // p.use_rho0 makes both rho0, which is what the asymmetry should be.
        const Real rho_ref = p.use_rho0 ? rho0 : Real(1000.0);
        w_fr += gamma_k * Hz * (rho_ref + rho) * orhoi;
        temp += gamma_k * (Frazil::Lhat + temp * (Real(1.0) - Frazil::r));
    } else if ((w_fr > Real(0.0)) && (temp > t_fr)) {
        Real gamma_k = Frazil::pfac * (t_fr - temp)
            / (Frazil::Lhat + temp * (Real(1.0) - Frazil::r)
               + Real(0.0543) * salt);
        // ...and rho0 here, which is the inconsistency.
        const Real delta_wfr = gamma_k * Hz * (rho0 + rho) * orhoi;
        if ((w_fr + delta_wfr) > Real(0.0)) {
            w_fr += delta_wfr;
        } else {
            gamma_k = -w_fr * p.IceRho / (Hz * (rho0 + rho));
            w_fr = Real(0.0);
        }
        temp += gamma_k * (Frazil::Lhat + temp * (Real(1.0) - Frazil::r));
    }
}

} // namespace REMORASeaIce

#endif
