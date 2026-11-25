#include <REMORA.H>
#include <REMORA_SeaIce.H>
#include <REMORA_IceEVP.H>
#include <REMORA_IceElastic.H>

using namespace amrex;

/**
 * \brief The elastic-viscous-plastic ice dynamics subcycle.
 *
 * A port of the ICE_MOMENTUM/ICE_EVP block of ROMS seaice.F:
 *
 *     liuol = liunw ; liunw = 3 - liuol
 *     dtevp = dtice / nEVP
 *     DO my_iEVP = 1, nEVP
 *        lieol = lienw ; lienw = 3 - lieol
 *        ice_evp      ! viscosities from the current strain
 *        ice_evp_sig  ! step the internal stress
 *        ice_elastic  ! step the momentum
 *     END DO
 *
 * ROMS indexes the two ELASTIC time levels with lieol/lienw into the same
 * two slots of Si that liold/linew and liuol/liunw also index, for different
 * components. Here each alternating pair gets its own storage and the loop
 * ping-pongs between two MultiFabs, which is the same arithmetic without
 * three index pairs aliasing one array.
 *
 * At the last subcycle the elastic velocity becomes the ice velocity, which
 * is what ice_vbc and ice_spdiw read on the next step.
 *
 * The lateral conditions are applied inside every subcycle, on the subcycle's
 * own old and new levels: to the three stresses after ice_evp_sig
 * (ice_evp_sig.F:181-205) and to the two velocities after ice_elastic
 * (ice_elastic.F:727-737). A radiation side therefore radiates at the elastic
 * step, and the momentum reads the stress the condition set.
 */
void
REMORA::ice_dynamics (int lev, Real dtice)
{
    BL_PROFILE("REMORA::ice_dynamics()");
    if (!sea_ice_params.do_momentum) { return; }

    const auto& p = sea_ice_params;
    const int nEVP = p.nEVP;
    const Real dtevp = dtice / static_cast<Real>(nEVP);
    const Real l_rho0 = solverChoice.rho0;
    const Real l_g    = solverChoice.g;
    const bool l_zmax = p.zmax_clamp;

    const BoxArray&            ba = vec_ice_f[lev]->boxArray();
    const DistributionMapping& dm = vec_ice_f[lev]->DistributionMap();
    const IntVect              ng = vec_ice_f[lev]->nGrowVect();

    // The elastic state, two levels each. ROMS's lieol/lienw.
    MultiFab sig_old(ba, dm, REMORASeaIce::IceSig::NumComps, ng);
    MultiFab sig_new(ba, dm, REMORASeaIce::IceSig::NumComps, ng);
    MultiFab uevp_old(convert(ba,IntVect(1,0,0)), dm, 1, ng);
    MultiFab uevp_new(convert(ba,IntVect(1,0,0)), dm, 1, ng);
    MultiFab vevp_old(convert(ba,IntVect(0,1,0)), dm, 1, ng);
    MultiFab vevp_new(convert(ba,IntVect(0,1,0)), dm, 1, ng);
    // The viscosities and strength, recomputed every subcycle.
    MultiFab visc(ba, dm, 3, ng);   // 0 = Pice, 1 = Bvis, 2 = Svis
    sig_old.setVal(Real(0.0)); sig_new.setVal(Real(0.0)); visc.setVal(Real(0.0));

    // ice_evp.F reloads the elastic velocity from the ice velocity on the
    // FIRST subcycle; everything after that carries over.
    MultiFab::Copy(uevp_old, *vec_ice_u[lev], 0,0,1,ng);
    MultiFab::Copy(vevp_old, *vec_ice_v[lev], 0,0,1,ng);
    MultiFab::Copy(uevp_new, *vec_ice_u[lev], 0,0,1,ng);
    MultiFab::Copy(vevp_new, *vec_ice_v[lev], 0,0,1,ng);
    // The stress carries over between ice steps, so it starts from Si.
    // sig_bc_var is the row of the boundary-condition table for each.
    int sig_bc_var[REMORASeaIce::IceSig::NumComps];
    for (int c = 0; c < REMORASeaIce::IceSig::NumComps; ++c) {
        const int src = (c == REMORASeaIce::IceSig::xx) ? REMORASeaIce::IceS::ISxx
                      : (c == REMORASeaIce::IceSig::yy) ? REMORASeaIce::IceS::ISyy
                                                        : REMORASeaIce::IceS::ISxy;
        MultiFab::Copy(sig_old, *vec_ice_s_new[lev], src, c, 1, ng);
        sig_bc_var[c] = REMORASeaIce::SeaIceParameters::ice_bc_var_for_comp(src);
    }

    for (int sub = 1; sub <= nEVP; ++sub)
    {
        // ---- ice_evp: viscosities from the current strain ----
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(visc, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Array4<Real      > const& vv = visc.array(mfi);
            Array4<Real const> const& Si = vec_ice_s_new[lev]->const_array(mfi);
            Array4<Real const> const& uo = uevp_old.const_array(mfi);
            Array4<Real const> const& vo = vevp_old.const_array(mfi);
            Array4<Real const> const& pm = vec_pm[lev]->const_array(mfi);
            Array4<Real const> const& pn = vec_pn[lev]->const_array(mfi);
            Array4<Real const> const& Fi = vec_ice_f[lev]->const_array(mfi);
            // Interior points, which are all the stress below reads.
            Box bx = mfi.tilebox();
            bx.setSmall(2,0); bx.setBig(2,0);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                const auto e = REMORASeaIce::ice_strain_rates(
                    uo(i,j,0), uo(i+1,j,0), uo(i,j+1,0), uo(i+1,j+1,0),
                    uo(i,j-1,0), uo(i+1,j-1,0),
                    vo(i,j,0), vo(i,j+1,0), vo(i+1,j,0), vo(i+1,j+1,0),
                    vo(i-1,j,0), vo(i-1,j+1,0),
                    pm(i,j,0), pn(i,j,0));
                Real P, B, S;
                REMORASeaIce::ice_evp_viscosity(
                    e, Si(i,j,0,REMORASeaIce::IceS::Hice),
                    Si(i,j,0,REMORASeaIce::IceS::Aice), p, l_zmax, P, B, S,
                    Fi(i,j,0,REMORASeaIce::IceF::Pgrd));
                vv(i,j,0,0) = P; vv(i,j,0,1) = B; vv(i,j,0,2) = S;
            });
        }

        // ---- ice_evp_sig: step the internal stress ----
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(sig_new, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Array4<Real      > const& sn = sig_new.array(mfi);
            Array4<Real const> const& so = sig_old.const_array(mfi);
            Array4<Real const> const& vv = visc.const_array(mfi);
            Array4<Real const> const& Si = vec_ice_s_new[lev]->const_array(mfi);
            Array4<Real const> const& uo = uevp_old.const_array(mfi);
            Array4<Real const> const& vo = vevp_old.const_array(mfi);
            Array4<Real const> const& pm = vec_pm[lev]->const_array(mfi);
            Array4<Real const> const& pn = vec_pn[lev]->const_array(mfi);
            // Interior points only, as ice_evp_sig.F: the ghost points of the
            // new stress come from the exchange and the conditions below.
            Box bx = mfi.tilebox();
            bx.setSmall(2,0); bx.setBig(2,0);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                const auto e = REMORASeaIce::ice_strain_rates(
                    uo(i,j,0), uo(i+1,j,0), uo(i,j+1,0), uo(i+1,j+1,0),
                    uo(i,j-1,0), uo(i+1,j-1,0),
                    vo(i,j,0), vo(i,j+1,0), vo(i+1,j,0), vo(i+1,j+1,0),
                    vo(i-1,j,0), vo(i-1,j+1,0),
                    pm(i,j,0), pn(i,j,0));
                Real xx, yy, xy;
                REMORASeaIce::ice_evp_stress(
                    e, Si(i,j,0,REMORASeaIce::IceS::Hice),
                    vv(i,j,0,0), vv(i,j,0,1), vv(i,j,0,2),
                    so(i,j,0,REMORASeaIce::IceSig::xx),
                    so(i,j,0,REMORASeaIce::IceSig::yy),
                    so(i,j,0,REMORASeaIce::IceSig::xy),
                    pm(i,j,0), dtevp, p, xx, yy, xy);
                sn(i,j,0,REMORASeaIce::IceSig::xx) = xx;
                sn(i,j,0,REMORASeaIce::IceSig::yy) = yy;
                sn(i,j,0,REMORASeaIce::IceSig::xy) = xy;
            });
        }
        sig_new.FillBoundary(geom[lev].periodicity());
        // ice_evp_sig.F:181-205 applies ice_bc2d to the three stresses in
        // every subcycle, with the subcycle's own old and new levels (lieol,
        // lienw) and the elastic velocity at lienw, which at this point still
        // holds the velocity of two subcycles ago. The momentum below reads
        // the boundary row one point beyond each box, hence the 1.
        for (const int c : {REMORASeaIce::IceSig::xx, REMORASeaIce::IceSig::xy,
                            REMORASeaIce::IceSig::yy}) {
            ice_bc2d(lev, sig_bc_var[c], sig_new, sig_old, c, uevp_new, vevp_new, 1);
        }

        // ---- ice_elastic: step the momentum ----
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(uevp_new, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Array4<Real      > const& un = uevp_new.array(mfi);
            Array4<Real      > const& vn = vevp_new.array(mfi);
            Array4<Real const> const& uo = uevp_old.const_array(mfi);
            Array4<Real const> const& vo = vevp_old.const_array(mfi);
            Array4<Real const> const& sn = sig_new.const_array(mfi);
            Array4<Real const> const& Si = vec_ice_s_new[lev]->const_array(mfi);
            Array4<Real const> const& Fi = vec_ice_f[lev]->const_array(mfi);
            Array4<Real const> const& ua = vec_ice_uavg[lev]->const_array(mfi);
            Array4<Real const> const& va = vec_ice_vavg[lev]->const_array(mfi);
            Array4<Real const> const& au = vec_ice_aius[lev]->const_array(mfi);
            Array4<Real const> const& av = vec_ice_aivs[lev]->const_array(mfi);
            Array4<Real const> const& mu = vec_msku[lev]->const_array(mfi);
            Array4<Real const> const& mv = vec_mskv[lev]->const_array(mfi);
            Array4<Real const> const& pm = vec_pm[lev]->const_array(mfi);
            Array4<Real const> const& pn = vec_pn[lev]->const_array(mfi);
            Array4<Real const> const& fc = vec_fcor[lev]->const_array(mfi);
            Box ubx = mfi.nodaltilebox(0); ubx.setSmall(2,0); ubx.setBig(2,0);
            Box vbx = mfi.nodaltilebox(1); vbx.setSmall(2,0); vbx.setBig(2,0);
            ParallelFor(ubx, [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                un(i,j,0) = REMORASeaIce::ice_elastic_u(
                    i, j, Si, sn, Fi, uo, vo, ua, au, mu, mv, pm, pn,
                    fc, dtevp, l_rho0, l_g, p);
            });
            ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                vn(i,j,0) = REMORASeaIce::ice_elastic_v(
                    i, j, Si, sn, Fi, uo, vo, va, av, mu, mv, pm, pn,
                    fc, dtevp, l_rho0, l_g, p);
            });
        }
        // ice_elastic.F applies the velocity conditions inside the subcycle, on
        // the subcycle's own old and new levels. The step passed is only the
        // nudging's, ROMS's `tau*dt(ng)`; the radiation has no time step in it.
        // A face that two boxes share takes one value before the conditions read
        // the new velocity and again for what they write.
        uevp_new.OverrideSync(geom[lev].periodicity());
        vevp_new.OverrideSync(geom[lev].periodicity());
        ice_velbc(lev, uevp_new, vevp_new, uevp_old, vevp_old, dtice);
        uevp_new.OverrideSync(geom[lev].periodicity());
        vevp_new.OverrideSync(geom[lev].periodicity());

        uevp_new.FillBoundary(geom[lev].periodicity());
        vevp_new.FillBoundary(geom[lev].periodicity());

        // Ping-pong the elastic levels, as ROMS flips lieol/lienw.
        std::swap(sig_old,  sig_new);
        std::swap(uevp_old, uevp_new);
        std::swap(vevp_old, vevp_new);
    }

    // The last subcycle's elastic velocity becomes the ice velocity, and the
    // stress is carried back into Si for the next ice step. After the final
    // swap the newest values are in the *_old handles.
    // Keep the velocity this step started from: ROMS's ice_spdiw and ice_vbc
    // index the ice velocity with liuol BEFORE seaice.F swaps liuol and liunw,
    // so on the next step they read this one, not the one written below.
    MultiFab::Copy(*vec_ice_u_prev[lev], *vec_ice_u[lev], 0,0,1,ng);
    MultiFab::Copy(*vec_ice_v_prev[lev], *vec_ice_v[lev], 0,0,1,ng);
    MultiFab::Copy(*vec_ice_u[lev], uevp_old, 0,0,1,ng);
    MultiFab::Copy(*vec_ice_v[lev], vevp_old, 0,0,1,ng);
    // The elastic velocity offered to the history file is the last subcycle's
    MultiFab::Copy(*vec_ice_uevp[lev], uevp_old, 0,0,1,ng);
    MultiFab::Copy(*vec_ice_vevp[lev], vevp_old, 0,0,1,ng);
    for (int c = 0; c < REMORASeaIce::IceSig::NumComps; ++c) {
        const int dst = (c == REMORASeaIce::IceSig::xx) ? REMORASeaIce::IceS::ISxx
                      : (c == REMORASeaIce::IceSig::yy) ? REMORASeaIce::IceS::ISyy
                                                        : REMORASeaIce::IceS::ISxy;
        // With its ghost points: they are the old level of the next ice
        // step's first stress condition.
        MultiFab::Copy(*vec_ice_s_new[lev], sig_old, c, dst, 1, ng);
    }
    vec_ice_u[lev]->FillBoundary(geom[lev].periodicity());
    vec_ice_v[lev]->FillBoundary(geom[lev].periodicity());
    vec_ice_s_new[lev]->FillBoundary(geom[lev].periodicity());
}
