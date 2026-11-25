#include <REMORA.H>
#include <REMORA_SeaIce.H>
#ifdef REMORA_USE_NETCDF
#include <REMORA_NCFile.H>
#endif

using namespace amrex;

/**
 * \file
 * \brief Time-varying ice open-boundary data: ROMS's BOUNDARY%iceG_* plus
 *        set_ngfld, on REMORA's own boundary time-series machinery.
 *
 * ROMS reads each boundary record with get_ngfld and interpolates it to the model time with
 * set_ngfld. Here NCTimeSeriesBoundary does both, as it does for the ocean.
 *
 * The records have ROMS's names, with the side appended as in varinfo.yaml: `Aice_west`,
 * `ice_thickness_east`, `Uice_north` and so on. varinfo.yaml has such entries for only seven
 * of the ice variables, and ROMS's `clamped` ice temperature reads an unassigned field index;
 * here every variable whose condition can read a record has one, `Tice_west` for example, and
 * takes the constant `remora.ice.bc_value.<var>` where the file holds none.
 *
 * Differences from ROMS: the interpolation weights are not rounded to the millisecond, and a
 * file needs at least two time records. A constant boundary value is `remora.ice.bc_value`.
 */

void
REMORA::init_ice_bdry_from_netcdf (int lev)
{
    BL_PROFILE("REMORA::init_ice_bdry_from_netcdf()");
    const auto& p = sea_ice_params;

    if (!REMORASeaIce::has_sea_ice(sea_ice_model)) { return; }

    int nvar = 0;
    const auto* vars = REMORASeaIce::SeaIceParameters::ice_bc_var_names(nvar);

    // Which sides of which variables ROMS would `acquire`.
    bool any_needs_data = false;
    for (int v = 0; v < nvar && !any_needs_data; ++v) {
        for (int sd = 0; sd < REMORASeaIce::IceSide::count; ++sd) {
            if (p.ice_bc_needs_data(v, sd)) { any_needs_data = true; break; }
        }
    }

    const bool have_file = !nc_ice_bdry_file.empty() && !nc_ice_bdry_file[0].empty();

#ifndef REMORA_USE_NETCDF
    amrex::ignore_unused(lev, vars, any_needs_data, have_file);
#else
    ice_boundary_series[lev].resize(nvar);

    if (!have_file) {
        // With every side `closed`, the default, no value is needed
        if (any_needs_data && lev == 0) {
            Print() << "Sea ice: open boundaries ask for data and "
                    << "remora.ice.nc_bdry_file is not set, so every "
                    << "boundary value comes from remora.ice.bc_value.*\n";
        }
        return;
    }

    int rx = 1, ry = 1;
    for (int k = lev - 1; k >= 0; --k) { rx *= ref_ratio[k][0]; ry *= ref_ratio[k][1]; }

    for (int v = 0; v < nvar; ++v) {
        GpuArray<bool, AMREX_SPACEDIM*2> need{};
        for (int ori = 0; ori < AMREX_SPACEDIM*2; ++ori) { need[ori] = false; }
        need[Orientation(Direction::x, Orientation::low)]  = p.ice_bc_needs_data(v, REMORASeaIce::IceSide::west);
        need[Orientation(Direction::x, Orientation::high)] = p.ice_bc_needs_data(v, REMORASeaIce::IceSide::east);
        need[Orientation(Direction::y, Orientation::low)]  = p.ice_bc_needs_data(v, REMORASeaIce::IceSide::south);
        need[Orientation(Direction::y, Orientation::high)] = p.ice_bc_needs_data(v, REMORASeaIce::IceSide::north);

        bool any_side = false;
        for (int ori = 0; ori < AMREX_SPACEDIM*2; ++ori) { any_side = any_side || need[ori]; }
        if (!any_side) { continue; }

        // A variable whose records are not all in the file takes the constant instead. The
        // choice is printed.
        Vector<std::string> want;
        if (need[Orientation(Direction::x, Orientation::low)])  { want.push_back(std::string(vars[v]) + "_west"); }
        if (need[Orientation(Direction::x, Orientation::high)]) { want.push_back(std::string(vars[v]) + "_east"); }
        if (need[Orientation(Direction::y, Orientation::low)])  { want.push_back(std::string(vars[v]) + "_south"); }
        if (need[Orientation(Direction::y, Orientation::high)]) { want.push_back(std::string(vars[v]) + "_north"); }

        bool in_file = true;
        for (const auto& fn : nc_ice_bdry_file) {
            if (!QueryNetCDFHasVars(fn, want)) { in_file = false; break; }
        }

        std::string msg = "Sea ice: boundary values for ";
        msg += vars[v];
        if (!in_file) {
            msg += " come from remora.ice.bc_value (the file has no";
            for (const auto& w : want) { msg += " " + w; }
            msg += ")\n";
            if (lev == 0) { Print() << msg; }
            continue;
        }

        const std::string tname = sea_ice_bc_time_name[v].empty()
            ? std::string(REMORASeaIce::SeaIceParameters::ice_bc_time_name(v))
            : sea_ice_bc_time_name[v];

        const int stag = REMORASeaIce::SeaIceParameters::ice_bc_stagger(v);
        const IntVect itype = (stag == 1) ? IntVect(1,0,0)
                            : (stag == 2) ? IntVect(0,1,0)
                                          : IntVect(0,0,0);

        msg += " come from the file, time axis '" + tname + "', records";
        for (const auto& w : want) { msg += " " + w; }
        msg += "\n";
        if (lev == 0) { Print() << msg; }

        ice_boundary_series[lev][v].reset(new NCTimeSeriesBoundary(
            lev, geom, nc_ice_bdry_file, std::string(vars[v]), tname,
            itype, need, /*is2d=*/true, rx, ry, solverChoice.time_ref, verbose));
        ice_boundary_series[lev][v]->Initialize();
    }
#endif
}

void
REMORA::update_ice_bdry_data (int lev, Real time)
{
    if (!REMORASeaIce::has_sea_ice(sea_ice_model)) { return; }
#ifdef REMORA_USE_NETCDF
    if (lev >= ice_boundary_series.size()) { return; }
    for (auto& s : ice_boundary_series[lev]) {
        if (s) { s->update_interpolated_to_time(model_time(time)); }
    }
#else
    amrex::ignore_unused(lev, time);
#endif
}
