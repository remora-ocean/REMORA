#ifndef REMORA_ICEELASTIC_H_
#define REMORA_ICEELASTIC_H_

#include <algorithm>
#include <cmath>

#include <AMReX_REAL.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>

#include <REMORA_SeaIce.H>
#include <REMORA_IceEVP.H>

namespace REMORASeaIce {

/**
 * The ice momentum equations, ported from ROMS ice_elastic.F.
 *
 * The stencil depends on the mask: the ISxy average switches between a 2-point and a 4-point
 * form on the neighbouring umask/vmask.
 *
 * The update is semi-implicit in the ice-ocean drag: the drag appears
 * explicitly in `uforce` through Uavg and implicitly in `alfa` through the
 * new Uevp.
 *
 * The Coriolis term is added in u and subtracted in v.
 */

/** u-momentum. Returns the new Uevp at the u point (i,j). */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
ice_elastic_u (const int i, const int j,
               const amrex::Array4<amrex::Real const>& Si,    //!< rho state (IceS comps)
               const amrex::Array4<amrex::Real const>& sig,   //!< rho stress (IceSig comps)
               const amrex::Array4<amrex::Real const>& Fi,    //!< rho internal (IceF comps)
               const amrex::Array4<amrex::Real const>& Uevp,  //!< u points, current sub-level
               const amrex::Array4<amrex::Real const>& Vevp,  //!< v points, current sub-level
               const amrex::Array4<amrex::Real const>& Uavg,  //!< u points
               const amrex::Array4<amrex::Real const>& AIus,  //!< u points
               const amrex::Array4<amrex::Real const>& umask,
               const amrex::Array4<amrex::Real const>& vmask,
               const amrex::Array4<amrex::Real const>& pm,
               const amrex::Array4<amrex::Real const>& pn,
               const amrex::Array4<amrex::Real const>& fcor,
               const amrex::Real dtevp,
               const amrex::Real rho0,
               const amrex::Real g,
               const SeaIceParameters& p)
{
    using amrex::Real;
    const int H = IceS::Hice, A = IceS::Aice;
    const int XX = IceSig::xx, XY = IceSig::xy;

    // ROMS metrics.F: om_u = 2/(pm(i-1)+pm(i)) and on_u = 2/(pn(i-1)+pn(i))
    const Real pmu = Real(0.5) * (pm(i,j,0) + pm(i-1,j,0));
    const Real pnu = Real(0.5) * (pn(i,j,0) + pn(i-1,j,0));
    const Real om_u = Real(1.0) / pmu;
    const Real on_u = Real(1.0) / pnu;
    const Real area = om_u * on_u;
    const Real chux = Real(0.5) * (Fi(i-1,j,0,IceF::IOmt) + Fi(i,j,0,IceF::IOmt));
    Real masu = Real(0.5) * (Si(i,j,0,H) + Si(i-1,j,0,H));
    masu = amrex::max(masu, Real(0.1)) * p.IceRho;
    const Real auf = amrex::max(Real(0.5) * (Si(i-1,j,0,A) + Si(i,j,0,A)),
                                Real(0.01));

    const Real s1 = (sig(i,j,0,XX) - sig(i-1,j,0,XX)) * pmu;

    // The mask-dependent ISxy average: next to a closed face ROMS drops to
    // the 2-point form rather than reach across it.
    Real f1, f2;
    if ((umask(i,j  ,0) > Real(0.0)) && (umask(i,j+1,0) < Real(1.0))) {
        f1 = Real(0.5) * (sig(i,j,0,XY) + sig(i-1,j,0,XY));
    } else {
        f1 = Real(0.25) * (sig(i  ,j  ,0,XY) + sig(i  ,j+1,0,XY)
                         + sig(i-1,j+1,0,XY) + sig(i-1,j  ,0,XY));
    }
    if ((umask(i,j  ,0) > Real(0.0)) && (umask(i,j-1,0) < Real(1.0))) {
        f2 = Real(0.5) * (sig(i,j,0,XY) + sig(i-1,j,0,XY));
    } else {
        f2 = Real(0.25) * (sig(i  ,j  ,0,XY) + sig(i-1,j  ,0,XY)
                         + sig(i-1,j-1,0,XY) + sig(i  ,j-1,0,XY));
    }
    const Real s2 = (f1 - f2) * pnu;

    Real uforce = (s1 + s2) * area;
    uforce += p.IceRho * AIus(i,j,0) * area;
    uforce -= g * masu * on_u
            * (Fi(i,j,0,IceF::Hsse) - Fi(i-1,j,0,IceF::Hsse));

    // Coriolis, as a 2x2 block of mass fluxes normalised by the number of
    // open v faces.
    Real fakt = Real(0.0);
    const Real dsum = vmask(i-1,j  ,0) + vmask(i  ,j  ,0)
                    + vmask(i-1,j+1,0) + vmask(i  ,j+1,0);
    if (dsum > Real(0.0)) { fakt = Real(1.0) / dsum; }
    const Real mfv11 = Real(0.5) * (Si(i-1,j-1,0,H) * fcor(i-1,j-1,0)
                                  + Si(i-1,j  ,0,H) * fcor(i-1,j  ,0))
                     * Vevp(i-1,j,0) * vmask(i-1,j,0);
    const Real mfv21 = Real(0.5) * (Si(i,j-1,0,H) * fcor(i,j-1,0)
                                  + Si(i,j  ,0,H) * fcor(i,j  ,0))
                     * Vevp(i,j,0) * vmask(i,j,0);
    const Real mfv12 = Real(0.5) * (Si(i-1,j  ,0,H) * fcor(i-1,j  ,0)
                                  + Si(i-1,j+1,0,H) * fcor(i-1,j+1,0))
                     * Vevp(i-1,j+1,0) * vmask(i-1,j+1,0);
    const Real mfv22 = Real(0.5) * (Si(i,j  ,0,H) * fcor(i,j  ,0)
                                  + Si(i,j+1,0,H) * fcor(i,j+1,0))
                     * Vevp(i,j+1,0) * vmask(i,j+1,0);
    uforce += fakt * p.IceRho * area * (mfv11 + mfv21 + mfv12 + mfv22);

    uforce = uforce / area + auf * chux * rho0 * Uavg(i,j,0);
    const Real alfa = masu + dtevp * auf * rho0 * chux;
    return ((masu * Uevp(i,j,0) + dtevp * uforce) / alfa) * umask(i,j,0);
}

/** v-momentum. Returns the new Vevp at the v point (i,j). */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
ice_elastic_v (const int i, const int j,
               const amrex::Array4<amrex::Real const>& Si,
               const amrex::Array4<amrex::Real const>& sig,
               const amrex::Array4<amrex::Real const>& Fi,
               const amrex::Array4<amrex::Real const>& Uevp,
               const amrex::Array4<amrex::Real const>& Vevp,
               const amrex::Array4<amrex::Real const>& Vavg,
               const amrex::Array4<amrex::Real const>& AIvs,
               const amrex::Array4<amrex::Real const>& umask,
               const amrex::Array4<amrex::Real const>& vmask,
               const amrex::Array4<amrex::Real const>& pm,
               const amrex::Array4<amrex::Real const>& pn,
               const amrex::Array4<amrex::Real const>& fcor,
               const amrex::Real dtevp,
               const amrex::Real rho0,
               const amrex::Real g,
               const SeaIceParameters& p)
{
    using amrex::Real;
    const int H = IceS::Hice, A = IceS::Aice;
    const int YY = IceSig::yy, XY = IceSig::xy;

    const Real pmv = Real(0.5) * (pm(i,j,0) + pm(i,j-1,0));
    const Real pnv = Real(0.5) * (pn(i,j,0) + pn(i,j-1,0));
    const Real om_v = Real(1.0) / pmv;
    const Real on_v = Real(1.0) / pnv;
    const Real area = om_v * on_v;
    Real masv = Real(0.5) * (Si(i,j,0,H) + Si(i,j-1,0,H));
    masv = amrex::max(masv, Real(0.1)) * p.IceRho;
    const Real avf = amrex::max(Real(0.5) * (Si(i,j-1,0,A) + Si(i,j,0,A)),
                                Real(0.01));
    const Real chuy = Real(0.5) * (Fi(i,j-1,0,IceF::IOmt) + Fi(i,j,0,IceF::IOmt));

    const Real s1 = (sig(i,j,0,YY) - sig(i,j-1,0,YY)) * pnv;

    Real f1, f2;
    if ((vmask(i  ,j,0) > Real(0.0)) && (vmask(i+1,j,0) < Real(1.0))) {
        f1 = Real(0.5) * (sig(i,j,0,XY) + sig(i,j-1,0,XY));
    } else {
        f1 = Real(0.25) * (sig(i  ,j  ,0,XY) + sig(i+1,j  ,0,XY)
                         + sig(i+1,j-1,0,XY) + sig(i  ,j-1,0,XY));
    }
    if ((vmask(i  ,j,0) > Real(0.0)) && (vmask(i-1,j,0) < Real(1.0))) {
        f2 = Real(0.5) * (sig(i,j,0,XY) + sig(i,j-1,0,XY));
    } else {
        f2 = Real(0.25) * (sig(i  ,j  ,0,XY) + sig(i-1,j  ,0,XY)
                         + sig(i-1,j-1,0,XY) + sig(i  ,j-1,0,XY));
    }
    const Real s2 = (f1 - f2) * pmv;

    Real vforce = (s1 + s2) * area;
    vforce += p.IceRho * AIvs(i,j,0) * area;
    vforce -= g * masv * om_v
            * (Fi(i,j,0,IceF::Hsse) - Fi(i,j-1,0,IceF::Hsse));

    Real fakt = Real(0.0);
    const Real dsum = umask(i  ,j-1,0) + umask(i+1,j-1,0)
                    + umask(i  ,j  ,0) + umask(i+1,j  ,0);
    if (dsum > Real(0.0)) { fakt = Real(1.0) / dsum; }
    const Real mfu11 = Real(0.5) * (Si(i-1,j-1,0,H) * fcor(i-1,j-1,0)
                                  + Si(i  ,j-1,0,H) * fcor(i  ,j-1,0))
                     * Uevp(i,j-1,0) * umask(i,j-1,0);
    const Real mfu21 = Real(0.5) * (Si(i  ,j-1,0,H) * fcor(i  ,j-1,0)
                                  + Si(i+1,j-1,0,H) * fcor(i+1,j-1,0))
                     * Uevp(i+1,j-1,0) * umask(i+1,j-1,0);
    const Real mfu12 = Real(0.5) * (Si(i-1,j,0,H) * fcor(i-1,j,0)
                                  + Si(i  ,j,0,H) * fcor(i  ,j,0))
                     * Uevp(i,j,0) * umask(i,j,0);
    const Real mfu22 = Real(0.5) * (Si(i  ,j,0,H) * fcor(i  ,j,0)
                                  + Si(i+1,j,0,H) * fcor(i+1,j,0))
                     * Uevp(i+1,j,0) * umask(i+1,j,0);
    vforce -= fakt * p.IceRho * area * (mfu11 + mfu21 + mfu12 + mfu22);

    vforce = vforce / area + avf * chuy * rho0 * Vavg(i,j,0);
    const Real alfa = masv + dtevp * avf * rho0 * chuy;
    return ((masv * Vevp(i,j,0) + dtevp * vforce) / alfa) * vmask(i,j,0);
}

} // namespace REMORASeaIce

#endif
