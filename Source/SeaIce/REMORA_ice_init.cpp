#include <REMORA.H>
#include <REMORA_SeaIce.H>

#ifdef REMORA_USE_NETCDF
#include <REMORA_NCFile.H>
#endif

using namespace amrex;

/**
 * \brief Read the initial sea-ice state from the NetCDF initial file.
 *
 * Reads the fifteen Si fields ROMS get_state reads,
 * under the varinfo.yaml names, and only those the file holds; the rest keep
 * the ice-free state. The internal Fi fields are not read by either model:
 * they start from zero and are set by the first step.
 *
 * Four of the fifteen change no result, in ROMS or here, and are read so that a ROMS
 * initial or restart file loads unchanged: the enthalpy, the age thickness, and the elastic
 * velocities, all of which are recomputed before they are used.
 */
void
REMORA::init_ice_from_netcdf (int lev)
{
    BL_PROFILE("REMORA::init_ice_from_netcdf()");
#ifdef REMORA_USE_NETCDF
    if (!REMORASeaIce::has_sea_ice(sea_ice_model)) { return; }

    struct Field { const char* name; MultiFab* mf; int comp; const MultiFab* mask; };
    MultiFab* S = vec_ice_s_new[lev].get();
    const MultiFab* mr = vec_mskr[lev].get();
    const Vector<Field> fields = {
        {"Aice",               S, REMORASeaIce::IceS::Aice, mr},
        {"ice_thickness",      S, REMORASeaIce::IceS::Hice, mr},
        {"meltpond_thickness", S, REMORASeaIce::IceS::Hmel, mr},
        {"snow_thickness",     S, REMORASeaIce::IceS::Hsno, mr},
        {"ice_age",            S, REMORASeaIce::IceS::Iage, mr},
        {"ice_Sxx",            S, REMORASeaIce::IceS::ISxx, mr},
        {"ice_Sxy",            S, REMORASeaIce::IceS::ISxy, mr},
        {"ice_Syy",            S, REMORASeaIce::IceS::ISyy, mr},
        {"Tice",               S, REMORASeaIce::IceS::Tice, mr},
        {"ice_enthalpy",       S, REMORASeaIce::IceS::Enth, mr},
        {"ice_age_thickness",  S, REMORASeaIce::IceS::Hage, mr},
        {"Uice", vec_ice_u[lev].get(),    0, vec_msku[lev].get()},
        {"Vice", vec_ice_v[lev].get(),    0, vec_mskv[lev].get()},
        {"Uevp", vec_ice_uevp[lev].get(), 0, vec_msku[lev].get()},
        {"Vevp", vec_ice_vevp[lev].get(), 0, vec_mskv[lev].get()},
    };

    std::string found;
    for (const auto& f : fields) {
        bool in_file = true;
        for (int idx = 0; idx < num_boxes_at_level[lev]; idx++) {
            in_file = in_file && QueryNetCDFHasVars(nc_init_file[lev][idx], {f.name});
        }
        if (!in_file) { continue; }
        found += std::string(" ") + f.name;

        const NCRegion region = nc_local_region(*f.mf);
        for (int idx = 0; idx < num_boxes_at_level[lev]; idx++) {
            NCLocalFab nc_fab;
            Vector<NCLocalFab*> nc_fabs = {&nc_fab};
            Vector<std::string> nc_names = {f.name};
            Vector<enum NC_Data_Dims_Type> nc_dims = {NC_Data_Dims_Type::Time_SN_WE};
            BuildLocalFABsFromNetCDFFile<FArrayBox,Real>(boxes_at_level[lev][idx], region,
                                                         nc_init_file[lev][idx],
                                                         nc_names, nc_dims, nc_fabs);
            for (MFIter mfi(*f.mf, false); mfi.isValid(); ++mfi) {
                nc_fab.copy_to((*f.mf)[mfi], 0, f.comp, 1);
            }
        }
        // Land holds the file's fill value, in the ring outside the domain too; the ice
        // state is zero there.
        MultiFab::Multiply(*f.mf, *f.mask, 0, f.comp, 1, IntVect(1,1,0));
        f.mf->FillBoundary(f.comp, 1, geom[lev].periodicity());
    }

    // The velocity ice_spdiw and ice_vbc read on the first step is the initial
    // one. The second copy of the state starts from the initial state too,
    // unless it is advanced separately: ROMS reads the file into one copy
    // only, and without an advection the other is never filled from it.
    if (sea_ice_params.advects() || !sea_ice_params.no_advect_two_states) {
        MultiFab::Copy(*vec_ice_s_old[lev], *S, 0, 0, S->nComp(), S->nGrowVect());
    }
    MultiFab::Copy(*vec_ice_u_prev[lev], *vec_ice_u[lev], 0, 0, 1, vec_ice_u[lev]->nGrowVect());
    MultiFab::Copy(*vec_ice_v_prev[lev], *vec_ice_v[lev], 0, 0, 1, vec_ice_v[lev]->nGrowVect());

    // A restart reads the ice from the checkpoint instead
    if (!restart_chkfile.empty()) { return; }
    if (found.empty()) {
        Print() << "Sea ice: the initial file holds no ice state, so the ice starts ice-free\n";
    } else {
        Print() << "Sea ice: initial state read from the initial file:" << found << "\n";
    }
#else
    amrex::ignore_unused(lev);
#endif
}

/**
 * \brief Set the initial sea-ice state from the problem's analytic hook.
 *
 * The counterpart of init_ice_from_netcdf for remora.ic_type = analytic.
 * Problems that define no ice leave the state as init_stuff set it, which is
 * ice-free, so this is a no-op for every case that has no ice.
 */
void
REMORA::init_ice_analytic (int lev)
{
    BL_PROFILE("REMORA::init_ice_analytic()");
    if (!REMORASeaIce::has_sea_ice(sea_ice_model)) { return; }

    prob->init_analytic_seaice(lev, geom[lev], solverChoice, *this,
                               *vec_ice_s_new[lev], *vec_ice_u[lev], *vec_ice_v[lev]);

    vec_ice_s_new[lev]->FillBoundary(geom[lev].periodicity());
    vec_ice_u[lev]->FillBoundary(geom[lev].periodicity());
    vec_ice_v[lev]->FillBoundary(geom[lev].periodicity());

    // Same as after the file read: the second copy of the state and the
    // previous velocity start from the initial state.
    if (sea_ice_params.advects() || !sea_ice_params.no_advect_two_states) {
        MultiFab::Copy(*vec_ice_s_old[lev], *vec_ice_s_new[lev], 0, 0,
                       vec_ice_s_new[lev]->nComp(), vec_ice_s_new[lev]->nGrowVect());
    }
    MultiFab::Copy(*vec_ice_u_prev[lev], *vec_ice_u[lev], 0, 0, 1, vec_ice_u[lev]->nGrowVect());
    MultiFab::Copy(*vec_ice_v_prev[lev], *vec_ice_v[lev], 0, 0, 1, vec_ice_v[lev]->nGrowVect());
}
