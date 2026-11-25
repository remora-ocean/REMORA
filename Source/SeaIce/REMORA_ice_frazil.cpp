#include <REMORA.H>
#include <REMORA_SeaIce.H>
#include <REMORA_IceFrazil.H>

using namespace amrex;

/**
 * \brief Frazil ice: form ice in supercooled water and warm it back.
 *
 * A port of ROMS ice_frazil.F, a routine ROMS contains and never calls. It runs at the end
 * of the step when remora.ice.do_frazil is set. It walks each column from the bottom up,
 * converts supercooling into ice, warms that water toward the freezing point, and melts the
 * accumulated frazil again in any layer above freezing. The result is Fi(icW_fr), which the
 * next step's thermodynamics reads as an ice source.
 */
void
REMORA::ice_frazil (int lev, Real dtice)
{
    BL_PROFILE("REMORA::ice_frazil()");
    if (!REMORASeaIce::has_sea_ice(sea_ice_model)) { return; }
    if (!sea_ice_params.do_frazil) { return; }

    const auto& p = sea_ice_params;
    const int N = geom[lev].Domain().bigEnd(2);
    const Real l_rho0 = solverChoice.rho0;
    const bool l_pfr  = p.pressure_freezing;

    MultiFab& S_new = *cons_new[lev];

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*vec_ice_f[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real      > const& Fi   = vec_ice_f[lev]->array(mfi);
        Array4<Real      > const& cons = S_new.array(mfi);
        Array4<Real const> const& z_r  = vec_z_r[lev]->const_array(mfi);
        Array4<Real const> const& Hz   = vec_Hz[lev]->const_array(mfi);
        Array4<Real const> const& rho  = vec_rho_ice[lev]->const_array(mfi);
        Array4<Real const> const& mskr = vec_mskr[lev]->const_array(mfi);

        Box bx = mfi.tilebox(); bx.setSmall(2,0); bx.setBig(2,0);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            Fi(i,j,0,REMORASeaIce::IceF::W_fr) = Real(0.0);
            // Land cells hold uninitialised memory on the GPU, and ROMS
            // skips them with its own overLand test.
            if (mskr(i,j,0) < Real(1.0)) { return; }

            Real w_fr = Real(0.0);
            // Bottom up, as ROMS does: the melting branch can only draw on
            // frazil the layers below it already formed.
            for (int k = 0; k <= N; ++k) {
                Real t = cons(i,j,k,Temp_comp);
                REMORASeaIce::ice_frazil_level(
                    t, cons(i,j,k,Salt_comp), z_r(i,j,k), Hz(i,j,k),
                    rho(i,j,k), l_rho0, (k == N), p, l_pfr, w_fr);
                cons(i,j,k,Temp_comp) = t;
            }
            w_fr /= dtice;
            Fi(i,j,0,REMORASeaIce::IceF::W_fr) =
                amrex::max(w_fr, Real(0.0)) * mskr(i,j,0);
        });
    }
    vec_ice_f[lev]->FillBoundary(geom[lev].periodicity());
}
