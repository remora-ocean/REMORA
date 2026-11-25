#include <REMORA.H>
#include <REMORA_SeaIce.H>
#include <REMORA_IceAdvect.H>

using namespace amrex;

namespace {

/**
 * \brief MPDATA one ice tracer, in place.
 *
 * ROMS ice_smolar.h's ice_mpdata_tile, called once per field by
 * ice_advect.F. Two passes: donor-cell upwind, then an antidiffusive
 * correction from a pseudo-velocity.
 *
 * Between the passes the point past each non-periodic domain edge takes the
 * first interior value of the upwind result, as in ice_smolar.h:778-816. The
 * boundary conditions on the advected field are ice_bc2d's job.
 */
void
mpdata_one (MultiFab& mf_state, const int comp,
            const MultiFab& mf_u, const MultiFab& mf_v,
            const MultiFab& mf_pm, const MultiFab& mf_pn,
            const MultiFab& mf_mskr, const MultiFab& mf_msku,
            const MultiFab& mf_mskv,
            const Geometry& geo, const Real dtice,
            const REMORASeaIce::SeaIceParameters& p)
{
    BL_PROFILE("mpdata_one()");
    const BoxArray&            ba = mf_state.boxArray();
    const DistributionMapping& dm = mf_state.DistributionMap();
    const IntVect              ng = mf_state.nGrowVect();

    MultiFab aif(ba, dm, 1, ng);
    MultiFab aflxu(convert(ba,IntVect(1,0,0)), dm, 1, ng);
    MultiFab aflxv(convert(ba,IntVect(0,1,0)), dm, 1, ng);
    MultiFab FE(convert(ba,IntVect(0,1,0)), dm, 1, ng);
    MultiFab FX(convert(ba,IntVect(1,0,0)), dm, 1, ng);
    aif.setVal(Real(0.0));   aflxu.setVal(Real(0.0)); aflxv.setVal(Real(0.0));
    FE.setVal(Real(0.0));    FX.setVal(Real(0.0));

    // ---- pass 1: donor-cell upwind ----
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(aif, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real      > const& fu = aflxu.array(mfi);
        Array4<Real      > const& fv = aflxv.array(mfi);
        Array4<Real const> const& fld = mf_state.const_array(mfi, comp);
        Array4<Real const> const& u  = mf_u.const_array(mfi);
        Array4<Real const> const& v  = mf_v.const_array(mfi);
        Array4<Real const> const& pm = mf_pm.const_array(mfi);
        Array4<Real const> const& pn = mf_pn.const_array(mfi);
        // on_u = 2/(pn(i-1)+pn(i)), om_v = 2/(pm(j-1)+pm(j)), as ROMS
        // metrics.F defines them. Derived rather than stored, so pm cannot be
        // handed in where a face length was meant.
        Box ubx = mfi.grownnodaltilebox(0, IntVect(1,1,0));
        Box vbx = mfi.grownnodaltilebox(1, IntVect(1,1,0));
        ubx.setSmall(2,0); ubx.setBig(2,0);
        vbx.setSmall(2,0); vbx.setBig(2,0);
        ParallelFor(ubx, [=] AMREX_GPU_DEVICE (int i, int j, int) {
            const Real on_u = Real(2.0) / (pn(i-1,j,0) + pn(i,j,0));
            fu(i,j,0) = REMORASeaIce::mpdata_upwind_u(u(i,j,0), fld(i-1,j,0),
                                                      fld(i,j,0), on_u);
        });
        ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int) {
            const Real om_v = Real(2.0) / (pm(i,j-1,0) + pm(i,j,0));
            fv(i,j,0) = REMORASeaIce::mpdata_upwind_v(v(i,j,0), fld(i,j-1,0),
                                                      fld(i,j,0), om_v);
        });
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(aif, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real      > const& a  = aif.array(mfi);
        Array4<Real const> const& fld = mf_state.const_array(mfi, comp);
        Array4<Real const> const& fu = aflxu.const_array(mfi);
        Array4<Real const> const& fv = aflxv.const_array(mfi);
        Array4<Real const> const& pm = mf_pm.const_array(mfi);
        Array4<Real const> const& pn = mf_pn.const_array(mfi);
        Array4<Real const> const& mr = mf_mskr.const_array(mfi);
        Box bx = mfi.growntilebox(IntVect(1,1,0));
        bx.setSmall(2,0); bx.setBig(2,0);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int) {
            // ar = 1/omn, and omn = 1/(pm*pn).
            const Real ar = pm(i,j,0) * pn(i,j,0);
            a(i,j,0) = (fld(i,j,0)
                      - dtice * (fu(i+1,j,0) - fu(i,j,0)
                               + fv(i,j+1,0) - fv(i,j,0)) * ar) * mr(i,j,0);
        });
    }
    aif.FillBoundary(geo.periodicity());

    // ROMS ICE_UPWIND: the upwind result is the answer.
    if (p.advect_upwind) {
        MultiFab::Copy(mf_state, aif, 0, comp, 1, 0);
        mf_state.FillBoundary(comp, 1, IntVect(1,1,0), geo.periodicity());
        return;
    }

    // The antidiffusive pass reads aif one cell past the domain. At a non-periodic
    // edge ROMS sets that cell to the first interior value, and a corner to the
    // corner cell, whatever the boundary condition; the copy is not re-masked.
    {
        const Box& dom = geo.Domain();
        const int ilo = dom.smallEnd(0), ihi = dom.bigEnd(0);
        const int jlo = dom.smallEnd(1), jhi = dom.bigEnd(1);
        const bool perx = geo.isPeriodic(0);
        const bool pery = geo.isPeriodic(1);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(aif, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Array4<Real> const& a = aif.array(mfi);
            Box bx = mfi.growntilebox(IntVect(1,1,0));
            bx.setSmall(2,0); bx.setBig(2,0);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int) {
                const bool outx = !perx && (i < ilo || i > ihi);
                const bool outy = !pery && (j < jlo || j > jhi);
                if (outx || outy) {
                    const int ii = perx ? i : amrex::Clamp(i, ilo, ihi);
                    const int jj = pery ? j : amrex::Clamp(j, jlo, jhi);
                    a(i,j,0) = a(ii,jj,0);
                }
            });
        }
    }

    // ---- pass 2: the masked differences the pseudo-velocity needs ----
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(aif, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real      > const& fe = FE.array(mfi);
        Array4<Real      > const& fx = FX.array(mfi);
        Array4<Real const> const& a  = aif.const_array(mfi);
        Array4<Real const> const& mu = mf_msku.const_array(mfi);
        Array4<Real const> const& mv = mf_mskv.const_array(mfi);
        Box ubx = mfi.grownnodaltilebox(0, IntVect(1,1,0));
        Box vbx = mfi.grownnodaltilebox(1, IntVect(1,1,0));
        ubx.setSmall(2,0); ubx.setBig(2,0);
        vbx.setSmall(2,0); vbx.setBig(2,0);
        ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int) {
            fe(i,j,0) = Real(0.5) * mv(i,j,0) * (a(i,j,0) - a(i,j-1,0));
        });
        ParallelFor(ubx, [=] AMREX_GPU_DEVICE (int i, int j, int) {
            fx(i,j,0) = Real(0.5) * mu(i,j,0) * (a(i,j,0) - a(i-1,j,0));
        });
    }
    FE.FillBoundary(geo.periodicity());
    FX.FillBoundary(geo.periodicity());

    // ---- pass 3: the antidiffusive correction ----
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(aif, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real      > const& fu = aflxu.array(mfi);
        Array4<Real      > const& fv = aflxv.array(mfi);
        Array4<Real const> const& a  = aif.const_array(mfi);
        Array4<Real const> const& fe = FE.const_array(mfi);
        Array4<Real const> const& fx = FX.const_array(mfi);
        Array4<Real const> const& u  = mf_u.const_array(mfi);
        Array4<Real const> const& v  = mf_v.const_array(mfi);
        Array4<Real const> const& pm = mf_pm.const_array(mfi);
        Array4<Real const> const& pn = mf_pn.const_array(mfi);
        Box ubx = mfi.nodaltilebox(0); ubx.setSmall(2,0); ubx.setBig(2,0);
        Box vbx = mfi.nodaltilebox(1); vbx.setSmall(2,0); vbx.setBig(2,0);
        ParallelFor(ubx, [=] AMREX_GPU_DEVICE (int i, int j, int) {
            const Real on_u = Real(2.0) / (pn(i-1,j,0) + pn(i,j,0));
            fu(i,j,0) = REMORASeaIce::mpdata_anti_u(i, j, a, fe, u, v, pm, pn,
                                                    on_u, dtice);
        });
        ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int) {
            const Real om_v = Real(2.0) / (pm(i,j-1,0) + pm(i,j,0));
            fv(i,j,0) = REMORASeaIce::mpdata_anti_v(i, j, a, fx, u, v, pm, pn,
                                                    om_v, dtice);
        });
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(aif, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real      > const& fld = mf_state.array(mfi, comp);
        Array4<Real const> const& a  = aif.const_array(mfi);
        Array4<Real const> const& fu = aflxu.const_array(mfi);
        Array4<Real const> const& fv = aflxv.const_array(mfi);
        Array4<Real const> const& pm = mf_pm.const_array(mfi);
        Array4<Real const> const& pn = mf_pn.const_array(mfi);
        Array4<Real const> const& mr = mf_mskr.const_array(mfi);
        Box bx = mfi.tilebox(); bx.setSmall(2,0); bx.setBig(2,0);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int) {
            fld(i,j,0) = (a(i,j,0)
                        - dtice * pm(i,j,0) * pn(i,j,0)
                        * (fu(i+1,j,0) - fu(i,j,0)
                         + fv(i,j+1,0) - fv(i,j,0))) * mr(i,j,0);
        });
    }
    mf_state.FillBoundary(comp, 1, IntVect(1,1,0), geo.periodicity());
}

} // namespace

/**
 * \brief Advect the ice tracers and apply the limiter.
 *
 * A port of ROMS ice_advect.F, which runs MPDATA once per field in a fixed
 * order, and of ice_limit.F, which seaice.F calls straight afterwards.
 *
 * The enthalpy is the one field that is NEGATIVE, and MPDATA is
 * positive-definite -- its `rate` denominators are MAX(epsil, <sum of the
 * field>), which only guards a non-negative field. So ROMS negates the
 * enthalpy, advects it, and negates it back, which is reproduced here.
 * Tice is then recovered from the advected enthalpy rather than advected
 * itself, because it is an intensive quantity.
 */
void
REMORA::ice_advect (int lev, Real dtice)
{
    BL_PROFILE("REMORA::ice_advect()");
    if (!sea_ice_params.advects()) { return; }

    const auto& p = sea_ice_params;
    using namespace REMORASeaIce;

    MultiFab& S = *vec_ice_s_new[lev];
    MultiFab Hice_before(S.boxArray(), S.DistributionMap(), 1, S.nGrowVect());
    MultiFab::Copy(Hice_before, S, IceS::Hice, 0, 1, S.nGrowVect());

    // ice_enthalpy.F: Enth = Hice*Tice, zeroed where the ice is thinner than
    // min_hi. Computed here rather than carried, because Tice is recovered
    // from it after the advection.
    if (p.do_thermo) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(S, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real> const& Si = S.array(mfi);
        // Every ghost cell: the upwind fluxes one face outside the box read the
        // field two cells out, and Hice, Tice and Iage are current there while a
        // product formed last step is not.
        Box bx = mfi.growntilebox(S.nGrowVect());
        bx.setSmall(2,0); bx.setBig(2,0);
        const Real min_hi = p.min_hi;
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int) {
            Si(i,j,0,IceS::Enth) = Si(i,j,0,IceS::Hice) * Si(i,j,0,IceS::Tice);
            if (Si(i,j,0,IceS::Hice) <= min_hi) { Si(i,j,0,IceS::Enth) = Real(0.0); }
            // ROMS negates the enthalpy for MPDATA; see the note above.
            Si(i,j,0,IceS::Enth) = -Si(i,j,0,IceS::Enth);
            // The age is advected as the thickness it belongs to (ice_smolar.h):
            // Hage = Hice*Iage from the state before the advection, zero where
            // that ice is thinner than min_hi. Iage is recovered afterwards.
            Si(i,j,0,IceS::Hage) = Si(i,j,0,IceS::Hice) * Si(i,j,0,IceS::Iage);
            if (Si(i,j,0,IceS::Hice) <= min_hi) { Si(i,j,0,IceS::Hage) = Real(0.0); }
        });
    }
    }

    // ice_smolar.h advects the snow, the melt ponds, the enthalpy and the age only under
    // ICE_THERMO.
    Vector<int> advected = {IceS::Aice, IceS::Hice};
    if (p.do_thermo) {
        advected.insert(advected.end(), {IceS::Hsno, IceS::Hmel, IceS::Enth, IceS::Hage});
    }
    for (int c : advected) {
        mpdata_one(S, c, *vec_ice_u[lev], *vec_ice_v[lev],
                   *vec_pm[lev], *vec_pn[lev], *vec_mskr[lev],
                   *vec_msku[lev], *vec_mskv[lev], geom[lev], dtice, p);
        // ice_smolar.h applies the lateral condition to each field as soon as it has
        // advected it. The enthalpy is negated here, so its condition waits until it is
        // restored below. The age thickness has none of its own: the condition is on
        // the age it is converted back into.
        if (c != IceS::Enth && c != IceS::Hage) { ice_bc2d(lev, c); }
    }

    // Iage = Hage/Hice from the advected pair, with ROMS's 1e-6 floor.
    if (p.do_thermo) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(S, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real> const& Si = S.array(mfi);
        Box bx = mfi.growntilebox(IntVect(1,1,0));
        bx.setSmall(2,0); bx.setBig(2,0);
        const Real min_hi = p.min_hi;
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int) {
            Si(i,j,0,IceS::Iage) = Si(i,j,0,IceS::Hage)
                / amrex::max(Si(i,j,0,IceS::Hice), Real(1.0e-6));
            if (Si(i,j,0,IceS::Hice) <= min_hi) {
                Si(i,j,0,IceS::Hage) = Real(0.0);
                Si(i,j,0,IceS::Iage) = Real(0.0);
            }
        });
    }
    ice_bc2d(lev, IceS::Iage);
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(S, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real      > const& Si = S.array(mfi);
        Array4<Real      > const& Fi = vec_ice_f[lev]->array(mfi);
        Array4<Real const> const& h0 = Hice_before.const_array(mfi);
        Box bx = mfi.growntilebox(IntVect(1,1,0));
        bx.setSmall(2,0); bx.setBig(2,0);
        const Real min_ai = p.min_ai;
        const auto pp = p;
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int) {
            // ice_advect.F: the ice divergence rate, from the thickness change
            // the advection produced.
            Fi(i,j,0,IceF::Wdiv) = (Si(i,j,0,IceS::Hice) - h0(i,j,0)) / dtice;

            // Undo the negation and recover Tice. 1e-6 is ROMS's floor on the
            // thickness in that division.
            if (pp.do_thermo) {
                Si(i,j,0,IceS::Enth) = -Si(i,j,0,IceS::Enth);
                Si(i,j,0,IceS::Tice) = Si(i,j,0,IceS::Enth)
                    / amrex::max(Si(i,j,0,IceS::Hice), Real(1.0e-6));
                if (Si(i,j,0,IceS::Aice) <= min_ai) {
                    Si(i,j,0,IceS::Enth) = Real(0.0);
                    Si(i,j,0,IceS::Tice) = Real(-2.0);
                }
            }

            // ice_limit.F, which seaice.F calls right after the advection.
            Real A = Si(i,j,0,IceS::Aice), H = Si(i,j,0,IceS::Hice);
            Real Sn = Si(i,j,0,IceS::Hsno), Hm = Si(i,j,0,IceS::Hmel);
            Real T = Si(i,j,0,IceS::Tice),  Ag = Si(i,j,0,IceS::Iage);
            ice_limit_state(A, H, Sn, Hm, T, Ag, pp);
            Si(i,j,0,IceS::Aice) = A;  Si(i,j,0,IceS::Hice) = H;
            Si(i,j,0,IceS::Hsno) = Sn; Si(i,j,0,IceS::Hmel) = Hm;
            Si(i,j,0,IceS::Tice) = T;  Si(i,j,0,IceS::Iage) = Ag;
        });
    }

    S.FillBoundary(geom[lev].periodicity());

    // ice_limit.F closes with the lateral conditions. The enthalpy and Tice are ice_tibc's.
    for (int c : {IceS::Aice, IceS::Hice, IceS::Hsno, IceS::Hmel, IceS::Iage}) {
        ice_bc2d(lev, c);
    }
    ice_tibc(lev);

    vec_ice_f[lev]->FillBoundary(geom[lev].periodicity());
}
