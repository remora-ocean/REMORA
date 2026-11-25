#include <REMORA.H>
#include <REMORA_SeaIce.H>
#include <REMORA_IceStress.H>

using namespace amrex;

/**
 * \brief Set a face field on the faces of a closed or open domain side.
 *
 * Corresponds to the normal-face part of ROMS bc_u2d_tile and bc_v2d_tile
 * (bc_2d.F), which ice_spdiw.F and ice_vbc.F apply to the fields they compute:
 * on a side that is not periodic, the face on the domain boundary is zero
 * where the barotropic velocity normal to it is closed, and takes the value of
 * the next face inward otherwise.
 *
 * @param[in   ] lev level to operate on
 * @param[inout] mf  field at u points (dir = 0) or v points (dir = 1)
 * @param[in   ] dir direction normal to the faces
 */
void
REMORA::ice_boundary_face_bc (int lev, MultiFab& mf, int dir)
{
    if (geom[lev].isPeriodic(dir)) { return; }

    const Box& domain = geom[lev].Domain();
    const int bc_idx = (dir == 0) ? ubar_bc() : vbar_bc();

    for (int side = 0; side < 2; ++side)
    {
        const Orientation ori(dir, (side == 0) ? Orientation::low : Orientation::high);
        const REMORA_BC type = phys_bc_type[bc_idx][ori];
        const bool closed = (type == REMORA_BC::slip_wall || type == REMORA_BC::no_slip_wall);
        // The boundary face, and the offset from it to the next face inward
        const int face  = (side == 0) ? domain.smallEnd(dir) : domain.bigEnd(dir) + 1;
        const int inward = (side == 0) ? 1 : -1;

        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            Box bx = mfi.validbox();
            if (bx.smallEnd(dir) > face || bx.bigEnd(dir) < face) { continue; }
            bx.setSmall(dir, face); bx.setBig(dir, face);
            bx.setSmall(2,0); bx.setBig(2,0);
            Array4<Real> const& a = mf.array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                if (closed) {
                    a(i,j,0) = Real(0.0);
                } else if (dir == 0) {
                    a(i,j,0) = a(i+inward,j,0);
                } else {
                    a(i,j,0) = a(i,j+inward,0);
                }
            });
        }
    }
}

namespace {
//! The first level below the depth mlio, searching down from the surface (ROMS ice_spdiw.F)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
int
ice_spdiw_bottom (int i, int j, int N, Real mlio, Array4<Real const> const& z_r)
{
    for (int k = N; k >= 0; --k) {
        if (z_r(i,j,k) < mlio) { return k; }
    }
    return 0;
}
}

/**
 * \brief Ice-ocean velocity shear and mixed-layer velocity.
 *
 * A port of ROMS ice_spdiw.F. Produces three internal fields the ice-ocean
 * stress needs:
 *
 *   Fi(icIOvs)  the magnitude of (ocean - ice) velocity at rho points
 *   Fi(icUavg)  thickness-weighted mean u over the top ~10 m, at u points
 *   Fi(icVavg)  the same for v
 *
 * The sum runs from the surface down to the first level whose center is deeper than ROMS's
 * hard-wired 10 m (without KPP mixing), clamped to the bottom level.
 */
void
REMORA::ice_spdiw (int lev)
{
    BL_PROFILE("REMORA::ice_spdiw()");
    if (!REMORASeaIce::has_sea_ice(sea_ice_model)) { return; }

    const int N = geom[lev].Domain().bigEnd(2);
    constexpr Real mlio = Real(-10.0);   // ROMS ice_spdiw.F

    MultiFab mf_uw(convert(vec_ice_f[lev]->boxArray(), IntVect(1,0,0)),
                   vec_ice_f[lev]->DistributionMap(), 1, vec_ice_f[lev]->nGrowVect());
    MultiFab mf_vw(convert(vec_ice_f[lev]->boxArray(), IntVect(0,1,0)),
                   vec_ice_f[lev]->DistributionMap(), 1, vec_ice_f[lev]->nGrowVect());
    mf_uw.setVal(Real(0.0));
    mf_vw.setVal(Real(0.0));

    // uw and vw first: the shear below reads them at i and i+1.
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*vec_ice_f[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real      > const& uw   = mf_uw.array(mfi);
        Array4<Real      > const& vw   = mf_vw.array(mfi);
        // ROMS u(:,:,:,nrhs), the start-of-step velocity
        Array4<Real const> const& u    = xvel_old[lev]->const_array(mfi);
        Array4<Real const> const& v    = yvel_old[lev]->const_array(mfi);
        Array4<Real const> const& z_r  = vec_z_r[lev]->const_array(mfi);
        Array4<Real const> const& z_w  = vec_z_w[lev]->const_array(mfi);

        Box ubx = mfi.nodaltilebox(0); ubx.setSmall(2,0); ubx.setBig(2,0);
        Box vbx = mfi.nodaltilebox(1); vbx.setSmall(2,0); vbx.setBig(2,0);

        ParallelFor(ubx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            const int nb = static_cast<int>(amrex::Math::round(
                Real(0.5) * (ice_spdiw_bottom(i-1,j,N,mlio,z_r) + ice_spdiw_bottom(i,j,N,mlio,z_r))));
            const int nbotu = amrex::max(amrex::min(nb, N), 0);
            Real acc = Real(0.0), totml = Real(0.0);
            for (int k = N; k >= nbotu; --k) {
                const Real dml = Real(0.5) * ((z_w(i-1,j,k+1) - z_w(i-1,j,k))
                                            + (z_w(i  ,j,k+1) - z_w(i  ,j,k)));
                acc   += u(i,j,k) * dml;
                totml += dml;
            }
            uw(i,j,0) = (totml > Real(0.0)) ? (acc / totml) : Real(0.0);
        });

        ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            const int nb = static_cast<int>(amrex::Math::round(
                Real(0.5) * (ice_spdiw_bottom(i,j-1,N,mlio,z_r) + ice_spdiw_bottom(i,j,N,mlio,z_r))));
            const int nbotv = amrex::max(amrex::min(nb, N), 0);
            Real acc = Real(0.0), totml = Real(0.0);
            for (int k = N; k >= nbotv; --k) {
                const Real dml = Real(0.5) * ((z_w(i,j-1,k+1) - z_w(i,j-1,k))
                                            + (z_w(i,j  ,k+1) - z_w(i,j  ,k)));
                acc   += v(i,j,k) * dml;
                totml += dml;
            }
            vw(i,j,0) = (totml > Real(0.0)) ? (acc / totml) : Real(0.0);
        });
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*vec_ice_f[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real      > const& Fi   = vec_ice_f[lev]->array(mfi);
        Array4<Real const> const& uw   = mf_uw.const_array(mfi);
        Array4<Real const> const& vw   = mf_vw.const_array(mfi);
        Array4<Real      > const& Uavg = vec_ice_uavg[lev]->array(mfi);
        Array4<Real      > const& Vavg = vec_ice_vavg[lev]->array(mfi);
        // ROMS Si(liuol): the velocity before the last dynamics step
        Array4<Real const> const& Uice = vec_ice_u_prev[lev]->const_array(mfi);
        Array4<Real const> const& Vice = vec_ice_v_prev[lev]->const_array(mfi);

        Box bx = mfi.tilebox(); bx.setSmall(2,0); bx.setBig(2,0);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            const Real du = (uw(i  ,j,0) - Uice(i  ,j,0))
                          + (uw(i+1,j,0) - Uice(i+1,j,0));
            const Real dv = (vw(i,j  ,0) - Vice(i,j  ,0))
                          + (vw(i,j+1,0) - Vice(i,j+1,0));
            Fi(i,j,0,REMORASeaIce::IceF::IOvs) =
                Real(0.5) * std::sqrt(du*du + dv*dv);
        });
        Box ubx = mfi.nodaltilebox(0); ubx.setSmall(2,0); ubx.setBig(2,0);
        Box vbx = mfi.nodaltilebox(1); vbx.setSmall(2,0); vbx.setBig(2,0);
        ParallelFor(ubx, [=] AMREX_GPU_DEVICE (int i, int j, int) {
            Uavg(i,j,0) = uw(i,j,0);
        });
        ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int) {
            Vavg(i,j,0) = vw(i,j,0);
        });
    }
    ice_boundary_face_bc(lev, *vec_ice_uavg[lev], 0);
    ice_boundary_face_bc(lev, *vec_ice_vavg[lev], 1);
    vec_ice_f[lev]->FillBoundary(geom[lev].periodicity());
    vec_ice_uavg[lev]->FillBoundary(geom[lev].periodicity());
    vec_ice_vavg[lev]->FillBoundary(geom[lev].periodicity());
}

/**
 * \brief Air-ice and ice-ocean stress; REPLACES sustr/svstr under ice.
 *
 * A port of ROMS ice_vbc.F. This is the piece that makes ice felt by the
 * ocean's momentum: ROMS forms
 *
 *   sustr = aix*chux*(Uice - Uavg) + (1 - aix)*sustr_ao
 *
 * so at full cover the wind stress on the water is REPLACED by the ice-water
 * drag. Without it, a fully ice-covered cell is still driven by the raw wind,
 * which is the one physically wrong thing a thermodynamics-only ice model
 * would otherwise do -- which is why ROMS calls this unconditionally from
 * seaice.F, outside the ICE_MOMENTUM guard.
 *
 * `kappa` here is ROMS's local 0.4 and `z0ii` is 0.01 m -- note that is NOT
 * the 0.02 m of ice_mk.h, which uses the same name for a different value.
 */
void
REMORA::ice_vbc (int lev)
{
    BL_PROFILE("REMORA::ice_vbc()");
    if (!REMORASeaIce::has_sea_ice(sea_ice_model)) { return; }

    const auto& p = sea_ice_params;
    const int N = geom[lev].Domain().bigEnd(2);

    MultiFab mf_chuiw(vec_ice_f[lev]->boxArray(), vec_ice_f[lev]->DistributionMap(),
                      1, vec_ice_f[lev]->nGrowVect());
    MultiFab mf_chuai(vec_ice_f[lev]->boxArray(), vec_ice_f[lev]->DistributionMap(),
                      1, vec_ice_f[lev]->nGrowVect());
    MultiFab mf_utau (vec_ice_f[lev]->boxArray(), vec_ice_f[lev]->DistributionMap(),
                      1, vec_ice_f[lev]->nGrowVect());
    mf_chuiw.setVal(Real(0.0)); mf_chuai.setVal(Real(0.0)); mf_utau.setVal(Real(0.0));

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*vec_ice_f[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real      > const& chuiw = mf_chuiw.array(mfi);
        Array4<Real      > const& chuai = mf_chuai.array(mfi);
        Array4<Real      > const& utau  = mf_utau.array(mfi);
        Array4<Real const> const& Fi    = vec_ice_f[lev]->const_array(mfi);
        Array4<Real const> const& Si    = vec_ice_s_new[lev]->const_array(mfi);
        Array4<Real const> const& z_r   = vec_z_r[lev]->const_array(mfi);
        Array4<Real const> const& z_w   = vec_z_w[lev]->const_array(mfi);

        // ROMS runs this loop over Jstr-1..Jend, Istr-1..Iend so the u- and
        // v-point averages below have their low neighbour.
        Box bx = mfi.growntilebox(IntVect(1,1,0));
        bx.setSmall(2,0); bx.setBig(2,0);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            const Real dztop = z_w(i,j,N+1) - z_r(i,j,N);
            Real ut, chu, chua;
            // icIOmt carries the PREVIOUS step's chuiw, as in ROMS.
            REMORASeaIce::ice_vbc_coef(Fi(i,j,0,REMORASeaIce::IceF::IOvs),
                                       Si(i,j,0,REMORASeaIce::IceS::Hice),
                                       Si(i,j,0,REMORASeaIce::IceS::Aice),
                                       dztop,
                                       Fi(i,j,0,REMORASeaIce::IceF::IOmt),
                                       p, ut, chu, chua);
            utau(i,j,0) = ut; chuiw(i,j,0) = chu; chuai(i,j,0) = chua;
        });
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*vec_ice_f[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real      > const& Fi    = vec_ice_f[lev]->array(mfi);
        Array4<Real const> const& Si    = vec_ice_s_new[lev]->const_array(mfi);
        Array4<Real const> const& chuiw = mf_chuiw.const_array(mfi);
        Array4<Real const> const& chuai = mf_chuai.const_array(mfi);
        Array4<Real const> const& utau  = mf_utau.const_array(mfi);
        Array4<Real      > const& sustr = vec_sustr[lev]->array(mfi);
        Array4<Real      > const& svstr = vec_svstr[lev]->array(mfi);
        Array4<Real const> const& su_ao = vec_sustr_ao[lev]->const_array(mfi);
        Array4<Real const> const& sv_ao = vec_svstr_ao[lev]->const_array(mfi);
        Array4<Real const> const& su_ai = vec_sustr_ai[lev]->const_array(mfi);
        Array4<Real const> const& sv_ai = vec_svstr_ai[lev]->const_array(mfi);
        Array4<Real      > const& aius  = vec_ice_aius[lev]->array(mfi);
        Array4<Real      > const& aivs  = vec_ice_aivs[lev]->array(mfi);
        Array4<Real const> const& Uavg  = vec_ice_uavg[lev]->const_array(mfi);
        Array4<Real const> const& Vavg  = vec_ice_vavg[lev]->const_array(mfi);
        Array4<Real const> const& Zt    = vec_Zt_avg1[lev]->const_array(mfi);
        // ROMS Si(liuol): the velocity before the last dynamics step
        Array4<Real const> const& Uice  = vec_ice_u_prev[lev]->const_array(mfi);
        Array4<Real const> const& Vice  = vec_ice_v_prev[lev]->const_array(mfi);

        Box ubx = mfi.nodaltilebox(0); ubx.setSmall(2,0); ubx.setBig(2,0);
        Box vbx = mfi.nodaltilebox(1); vbx.setSmall(2,0); vbx.setBig(2,0);
        Box bx  = mfi.tilebox();       bx.setSmall(2,0);  bx.setBig(2,0);

        ParallelFor(ubx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            const Real aix  = Real(0.5) * (Si(i,j,0,REMORASeaIce::IceS::Aice)
                                         + Si(i-1,j,0,REMORASeaIce::IceS::Aice));
            const Real chux = Real(0.5) * (chuiw(i,j,0) + chuiw(i-1,j,0));
            const Real chuax= Real(0.5) * (chuai(i,j,0) + chuai(i-1,j,0));
            Real sw, soi;
            REMORASeaIce::ice_vbc_stress(aix, chux, chuax, Uice(i,j,0), Uavg(i,j,0),
                                         su_ao(i,j,0), su_ai(i,j,0), su_ai(i-1,j,0),
                                         p, sw, soi);
            aius(i,j,0) = soi;  sustr(i,j,0) = sw;
        });
        ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            const Real aiy  = Real(0.5) * (Si(i,j,0,REMORASeaIce::IceS::Aice)
                                         + Si(i,j-1,0,REMORASeaIce::IceS::Aice));
            const Real chuy = Real(0.5) * (chuiw(i,j,0) + chuiw(i,j-1,0));
            const Real chuay= Real(0.5) * (chuai(i,j,0) + chuai(i,j-1,0));
            Real sw, soi;
            REMORASeaIce::ice_vbc_stress(aiy, chuy, chuay, Vice(i,j,0), Vavg(i,j,0),
                                         sv_ao(i,j,0), sv_ai(i,j,0), sv_ai(i,j-1,0),
                                         p, sw, soi);
            aivs(i,j,0) = soi;  svstr(i,j,0) = sw;
        });
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            Fi(i,j,0,REMORASeaIce::IceF::Hsse) = Zt(i,j,0);
            Fi(i,j,0,REMORASeaIce::IceF::IOmt) = chuiw(i,j,0);
            Fi(i,j,0,REMORASeaIce::IceF::IOfv) = utau(i,j,0);
        });
    }

    // The thermodynamics averages the stress on the two faces of a cell into its
    // friction velocity, so the value on a boundary face reaches the first row
    // of cells.
    ice_boundary_face_bc(lev, *vec_ice_aius[lev], 0);
    ice_boundary_face_bc(lev, *vec_ice_aivs[lev], 1);
    ice_boundary_face_bc(lev, *vec_sustr[lev], 0);
    ice_boundary_face_bc(lev, *vec_svstr[lev], 1);
    vec_ice_f[lev]->FillBoundary(geom[lev].periodicity());
    vec_ice_aius[lev]->FillBoundary(geom[lev].periodicity());
    vec_ice_aivs[lev]->FillBoundary(geom[lev].periodicity());
    vec_sustr[lev]->FillBoundary(geom[lev].periodicity());
    vec_svstr[lev]->FillBoundary(geom[lev].periodicity());
}
