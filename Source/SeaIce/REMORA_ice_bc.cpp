#include <REMORA.H>
#include <REMORA_SeaIce.H>
#include <REMORA_IceBC.H>

using namespace amrex;

/**
 * \file
 * \brief Applying the ice lateral boundary conditions to the model's fields.
 *
 * The kernels live in REMORA_IceBC.H. This file decides which cells of which
 * MultiFab each condition touches, in a domain decomposition ROMS does not have.
 *
 * Three things follow from how ROMS applies them:
 *
 *  1. ROMS applies the four sides SEQUENTIALLY to the same array, west then
 *     east then south then north, and then the corners. A later side reads
 *     what an earlier one wrote. So the side loop is OUTSIDE the MFIter loop
 *     and the corners come last -- not four independent passes.
 *
 *  2. A corner cell reads the two ghost cells its own box just filled. Any
 *     box that owns a corner of the domain owns both edges that meet there,
 *     so this needs no communication; it does need the sides to have run.
 *
 *  3. ROMS's loop bounds are not all the interior. The tangential `closed`
 *     runs over the TOTAL range (Istr..IendT), one further than the others,
 *     so free slip fills the corner cell as well.
 *
 * As in ROMS, they are called from inside each stage: the thermodynamics, the limiter, each
 * advected field, the internal stress, and the velocity in the elastic subcycle.
 */

namespace {

//! Does this box own the named edge of the domain?
bool
owns_edge (const Box& vbx, const Box& dom, int side)
{
    switch (side) {
        case REMORASeaIce::IceSide::west:  return vbx.smallEnd(0) == dom.smallEnd(0);
        case REMORASeaIce::IceSide::east:  return vbx.bigEnd(0)   == dom.bigEnd(0);
        case REMORASeaIce::IceSide::south: return vbx.smallEnd(1) == dom.smallEnd(1);
        default:                           return vbx.bigEnd(1)   == dom.bigEnd(1);
    }
}

#ifdef REMORA_USE_NETCDF
//! The ice BC side index as an amrex::Orientation, for the boundary fabs.
amrex::Orientation
side_orientation (int side)
{
    using amrex::Orientation;
    using amrex::Direction;
    switch (side) {
        case REMORASeaIce::IceSide::west:  return Orientation(Direction::x, Orientation::low);
        case REMORASeaIce::IceSide::east:  return Orientation(Direction::x, Orientation::high);
        case REMORASeaIce::IceSide::south: return Orientation(Direction::y, Orientation::low);
        default:                           return Orientation(Direction::y, Orientation::high);
    }
}
#endif

} // namespace

/**
 * The interpolated boundary data for one side of one ice BC variable, or an
 * EMPTY Array4 when that variable's values come from the constant
 * `remora.ice.bc_value.<var>` instead.
 *
 * The fab boxes `NCTimeSeriesBoundary` builds already sit exactly on the
 * cells these conditions write -- the ghost column or row for a scalar, the
 * boundary FACE for a velocity normal to the side, the ghost row for one
 * tangential to it -- because the stagger is passed in as the index type. So
 * the caller indexes the fab at the same (i,j) it is writing, with no
 * offset bookkeeping of its own.
 */
amrex::Array4<const amrex::Real>
REMORA::ice_bdry_array (int lev, int bc_var, int side) const
{
    using amrex::Array4;
#ifdef REMORA_USE_NETCDF
    if (lev >= ice_boundary_series.size()) { return Array4<const amrex::Real>(); }
    if (bc_var < 0 || bc_var >= ice_boundary_series[lev].size()) { return Array4<const amrex::Real>(); }
    const auto& s = ice_boundary_series[lev][bc_var];
    if (!s) { return Array4<const amrex::Real>(); }
    return s->dat_interp(side_orientation(side)).const_array();
#else
    amrex::ignore_unused(lev, bc_var, side);
    return Array4<const amrex::Real>();
#endif
}

/**
 * \brief One ice scalar's lateral boundary condition: ROMS ice_bc2d.F.
 *
 * `comp` is the component of vec_ice_s_new, whose ghost cells are written.
 */
void
REMORA::ice_bc2d (int lev, int comp)
{
    if (!REMORASeaIce::has_sea_ice(sea_ice_model)) { return; }
    const int bc_var = REMORASeaIce::SeaIceParameters::ice_bc_var_for_comp(comp);
    if (bc_var < 0) { return; }
    // The elastic step reads the boundary row one point beyond each box (the thickness, in the
    // Coriolis term), so the condition is evaluated there too and needs the halo of the row
    // inside the boundary.
    MultiFab& S = *vec_ice_s_new[lev];
    S.FillBoundary(comp, 1, geom[lev].periodicity());
    ice_bc2d(lev, bc_var, S, *vec_ice_s_old[lev], comp, *vec_ice_u[lev], *vec_ice_v[lev], 1);
}

/**
 * \brief ROMS ice_bc2d.F on an explicit pair of time levels.
 *
 * ROMS calls ice_bc2d_tile with the pair of levels of whichever stage it is
 * in: liold/linew for the advected fields, and the elastic subcycle's own
 * lieol/lienw for the three internal stresses (ice_evp_sig.F:181-205), where
 * the velocity the `mixed` branch tests is the elastic one at lienw.
 *
 * @param[in   ] lev    level
 * @param[in   ] bc_var row of the bc_type table that carries the condition
 * @param[inout] S      new level, whose ghost cells are written
 * @param[in   ] S_old  old level, ROMS's `know`
 * @param[in   ] comp   component of S and S_old
 * @param[in   ] mf_ui  ice u-velocity the `mixed` branch tests
 * @param[in   ] mf_vi  ice v-velocity the `mixed` branch tests
 * @param[in   ] nhalo  ghost points to fill along the side beyond each box,
 *                      where the box meets another box or a periodic edge
 *
 * ROMS exchanges the whole boundary row after the condition
 * (ice_evp_sig.F:208-229), so a point of that row in a neighbouring tile's or
 * a periodic halo holds the condition's value. FillBoundary copies valid
 * points only and the boundary row is outside every box, so with `nhalo` the
 * condition is evaluated at those points directly. That reads the halo of the
 * row inside the boundary, which the caller must have filled.
 */
void
REMORA::ice_bc2d (int lev, int bc_var, MultiFab& S, const MultiFab& S_old,
                  int comp, const MultiFab& mf_ui, const MultiFab& mf_vi,
                  int nhalo)
{
    BL_PROFILE("REMORA::ice_bc2d()");

    const auto& p   = sea_ice_params;
    const Box& dom  = geom[lev].Domain();
    const auto& per = geom[lev].periodicity();
    const Real dtice = p.bc_dt;

    for (int side = 0; side < REMORASeaIce::IceSide::count; ++side)
    {
        const auto type = p.bc_type[bc_var][side];
        if (type == REMORASeaIce::IceBCType::periodic) { continue; }
        // A periodic direction is handled by FillBoundary, and ROMS's own
        // blocks are skipped on a periodic edge too.
        const bool is_x = REMORASeaIce::ice_side_is_x(side);
        if (per.isPeriodic(is_x ? 0 : 1)) { continue; }

        const Real bval   = p.bc_value[bc_var][side];
        const bool nudge  = p.bc_nudging[side];
        const Real tau_in = p.bc_tau_in, tau_out = p.bc_tau_out;
        const bool rad    = (type == REMORASeaIce::IceBCType::radiation);
        // Time-varying data for this side, if the file had it. An empty array
        // means the constant above is the value, which is what every `closed`
        // side and every file-less run uses.
        const auto bdat = ice_bdry_array(lev, bc_var, side);
        const bool from_file = (bdat.p != nullptr);

        for (MFIter mfi(S, false); mfi.isValid(); ++mfi)
        {
            const Box& vbx = mfi.validbox();
            if (!owns_edge(vbx, dom, side)) { continue; }

            Array4<Real      > const& f_new = S.array(mfi);
            Array4<Real const> const& f_old = S_old.const_array(mfi);
            Array4<Real const> const& ui    = mf_ui.const_array(mfi);
            Array4<Real const> const& vi    = mf_vi.const_array(mfi);
            Array4<Real const> const& mskr  = vec_mskr[lev]->const_array(mfi);
            Array4<Real const> const& msku  = vec_msku[lev]->const_array(mfi);
            Array4<Real const> const& mskv  = vec_mskv[lev]->const_array(mfi);

            // The boundary-adjacent interior cell, and the sweep along the
            // side clipped to this box's interior range.
            const int ib = (side == REMORASeaIce::IceSide::west)  ? dom.smallEnd(0)
                         : (side == REMORASeaIce::IceSide::east)  ? dom.bigEnd(0) : 0;
            const int jb = (side == REMORASeaIce::IceSide::south) ? dom.smallEnd(1)
                         : (side == REMORASeaIce::IceSide::north) ? dom.bigEnd(1) : 0;
            // Along the side: this box's part of it, and `nhalo` further
            // where the box does not end at a physical edge.
            const int  t     = is_x ? 1 : 0;
            const bool t_per = per.isPeriodic(t);
            const int  tlo   = dom.smallEnd(t), tlen = dom.length(t);
            int lo = std::max(vbx.smallEnd(t), dom.smallEnd(t));
            int hi = std::min(vbx.bigEnd(t), dom.bigEnd(t));
            if (t_per || vbx.smallEnd(t) > dom.smallEnd(t)) { lo -= nhalo; }
            if (t_per || vbx.bigEnd(t)   < dom.bigEnd(t))   { hi += nhalo; }
            Box line(IntVect(is_x ? ib : lo, is_x ? lo : jb, 0),
                     IntVect(is_x ? ib : hi, is_x ? hi : jb, 0));

            const int l_side = side, l_comp = comp;
            const auto l_type = type;
            ParallelFor(line, [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                const int d  = REMORASeaIce::ice_side_dir(l_side);
                const int ig = is_x ? (i - d) : i;
                const int jg = is_x ? j       : (j - d);
                // ROMS's BOUNDARY(ng)%ice_west(ifld,j): one value per
                // boundary point, at the cell this condition writes. A point
                // in a periodic halo takes the value of the point it images.
                Real bv = bval;
                if (from_file) {
                    int tb = is_x ? jg : ig;
                    if (t_per) { tb = tlo + ((tb - tlo) % tlen + tlen) % tlen; }
                    bv = is_x ? bdat(ig,tb,0) : bdat(tb,jg,0);
                }
                Real out;
                if (rad) {
                    out = REMORASeaIce::ice_bc_radiation(
                        i, j, l_side, f_old, f_new, msku, mskv, mskr,
                        nudge, bv, tau_in, tau_out, dtice, l_comp);
                } else {
                    // ROMS's `mixed` switches on the velocity INTO the
                    // domain at the boundary face: +u on the west, -u on the
                    // east, +v on the south, -v on the north.
                    const Real vel_in =
                        (l_side == REMORASeaIce::IceSide::west)  ?  ui(i,j,0)
                      : (l_side == REMORASeaIce::IceSide::east)  ? -ui(i+1,j,0)
                      : (l_side == REMORASeaIce::IceSide::south) ?  vi(i,j,0)
                      :                                            -vi(i,j+1,0);
                    out = REMORASeaIce::ice_bc_simple(
                        i, j, l_side, l_type, f_old, f_new, mskr,
                        bv, vel_in, l_comp);
                }
                f_new(ig,jg,0,l_comp) = out;
            });
        }
    }

    // The corners, after every side, on whichever box owns each one.
    if (!per.isPeriodic(0) && !per.isPeriodic(1)) {
        for (MFIter mfi(S, false); mfi.isValid(); ++mfi)
        {
            const Box& vbx = mfi.validbox();
            Array4<Real> const& f_new = S.array(mfi);
            const int cx[2] = {dom.smallEnd(0) - 1, dom.bigEnd(0) + 1};
            const int cy[2] = {dom.smallEnd(1) - 1, dom.bigEnd(1) + 1};
            const int dx[2] = {+1, -1};
            const int l_comp = comp;
            for (int sx = 0; sx < 2; ++sx) {
                for (int sy = 0; sy < 2; ++sy) {
                    if (!owns_edge(vbx, dom, sx ? REMORASeaIce::IceSide::east
                                                : REMORASeaIce::IceSide::west)) { continue; }
                    if (!owns_edge(vbx, dom, sy ? REMORASeaIce::IceSide::north
                                                : REMORASeaIce::IceSide::south)) { continue; }
                    const int ci = cx[sx], cj = cy[sy];
                    const int di = dx[sx], dj = dx[sy];
                    Box corner(IntVect(ci,cj,0), IntVect(ci,cj,0));
                    ParallelFor(corner, [=] AMREX_GPU_DEVICE (int i, int j, int)
                    {
                        f_new(i,j,0,l_comp) = REMORASeaIce::ice_bc_corner(
                            i, j, di, dj, Array4<Real const>(f_new), l_comp);
                    });
                }
            }
        }
    }
}

/**
 * \brief The ice temperature and enthalpy boundary condition: ROMS ice_tibc.F.
 *
 * Not ice_bc2d applied to Tice: the condition is carried on the enthalpy
 * hi*ti and the temperature recovered from it, and there is no radiation
 * branch. Writes both IceS::Tice and IceS::Enth.
 */
void
REMORA::ice_tibc (int lev)
{
    BL_PROFILE("REMORA::ice_tibc()");
    if (!REMORASeaIce::has_sea_ice(sea_ice_model)) { return; }

    const auto& p   = sea_ice_params;
    const Box& dom  = geom[lev].Domain();
    const auto& per = geom[lev].periodicity();
    const Real min_hi = p.min_hi;
    const auto fix = p.bc_fix;

    namespace IceS = REMORASeaIce::IceS;
    const int bc_var  = REMORASeaIce::SeaIceParameters::ice_bc_var_for_comp(IceS::Tice);
    const int bc_hice = REMORASeaIce::SeaIceParameters::ice_bc_var_for_comp(IceS::Hice);

    MultiFab& S     = *vec_ice_s_new[lev];
    MultiFab& S_old = *vec_ice_s_old[lev];

    for (int side = 0; side < REMORASeaIce::IceSide::count; ++side)
    {
        const auto type = p.bc_type[bc_var][side];
        if (type == REMORASeaIce::IceBCType::periodic) { continue; }
        const bool is_x = REMORASeaIce::ice_side_is_x(side);
        if (per.isPeriodic(is_x ? 0 : 1)) { continue; }

        // The clamped branch takes the PRODUCT of the boundary thickness and
        // the boundary temperature, so it needs both values -- and they come
        // from two different records, which may independently be in the file
        // or not. The thickness is keyed on ITS OWN side's condition, not on
        // the temperature's, because ROMS reads whatever is in
        // BOUNDARY%ice_west(isHice,j) regardless of what the thickness's own
        // LBC says; so the thickness series is consulted whenever it exists.
        const Real b_hice = p.bc_value[bc_hice][side];
        const Real b_tice = p.bc_value[bc_var][side];
        const auto bdat_t = ice_bdry_array(lev, bc_var,  side);
        const auto bdat_h = ice_bdry_array(lev, bc_hice, side);
        const bool t_from_file = (bdat_t.p != nullptr);
        const bool h_from_file = (bdat_h.p != nullptr);

        for (MFIter mfi(S, false); mfi.isValid(); ++mfi)
        {
            const Box& vbx = mfi.validbox();
            if (!owns_edge(vbx, dom, side)) { continue; }

            Array4<Real      > const& s_new = S.array(mfi);
            Array4<Real const> const& s_old = S_old.const_array(mfi);
            Array4<Real const> const& ui    = vec_ice_u[lev]->const_array(mfi);
            Array4<Real const> const& vi    = vec_ice_v[lev]->const_array(mfi);
            Array4<Real const> const& mskr  = vec_mskr[lev]->const_array(mfi);

            const int ib = (side == REMORASeaIce::IceSide::west)  ? dom.smallEnd(0)
                         : (side == REMORASeaIce::IceSide::east)  ? dom.bigEnd(0) : 0;
            const int jb = (side == REMORASeaIce::IceSide::south) ? dom.smallEnd(1)
                         : (side == REMORASeaIce::IceSide::north) ? dom.bigEnd(1) : 0;
            const int lo = is_x ? std::max(vbx.smallEnd(1), dom.smallEnd(1))
                                : std::max(vbx.smallEnd(0), dom.smallEnd(0));
            const int hi = is_x ? std::min(vbx.bigEnd(1), dom.bigEnd(1))
                                : std::min(vbx.bigEnd(0), dom.bigEnd(0));
            Box line(IntVect(is_x ? ib : lo, is_x ? lo : jb, 0),
                     IntVect(is_x ? ib : hi, is_x ? hi : jb, 0));

            const int l_side = side;
            const auto l_type = type;
            ParallelFor(line, [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                const int d  = REMORASeaIce::ice_side_dir(l_side);
                const int ig = is_x ? (i - d) : i;
                const int jg = is_x ? j       : (j - d);
                const Real vel_in =
                    (l_side == REMORASeaIce::IceSide::west)  ?  ui(i,j,0)
                  : (l_side == REMORASeaIce::IceSide::east)  ? -ui(i+1,j,0)
                  : (l_side == REMORASeaIce::IceSide::south) ?  vi(i,j,0)
                  :                                            -vi(i,j+1,0);
                // Pre-loaded with the ghost cell's current values: a `mixed`
                // INFLOW on the east or south leaves the temperature alone,
                // which is ROMS's behaviour (see ice_tibc_cell).
                Real e  = s_new(ig,jg,0,REMORASeaIce::IceS::Enth);
                Real tt = s_new(ig,jg,0,REMORASeaIce::IceS::Tice);
                const Real bh = h_from_file ? bdat_h(ig,jg,0) : b_hice;
                const Real bt = t_from_file ? bdat_t(ig,jg,0) : b_tice;
                REMORASeaIce::ice_tibc_cell(
                    i, j, l_side, l_type,
                    s_old, s_old, s_new, s_new, s_old, s_new, mskr,
                    bh, bt, vel_in, min_hi, fix, e, tt,
                    REMORASeaIce::IceS::Hice, REMORASeaIce::IceS::Tice,
                    REMORASeaIce::IceS::Enth);
                s_new(ig,jg,0,REMORASeaIce::IceS::Enth) = e;
                s_new(ig,jg,0,REMORASeaIce::IceS::Tice) = tt;
            });
        }
    }

    if (!per.isPeriodic(0) && !per.isPeriodic(1)) {
        for (MFIter mfi(S, false); mfi.isValid(); ++mfi)
        {
            const Box& vbx = mfi.validbox();
            Array4<Real> const& s_new = S.array(mfi);
            Array4<Real const> const& mskr = vec_mskr[lev]->const_array(mfi);
            const int cx[2] = {dom.smallEnd(0) - 1, dom.bigEnd(0) + 1};
            const int cy[2] = {dom.smallEnd(1) - 1, dom.bigEnd(1) + 1};
            const int dx[2] = {+1, -1};
            for (int sx = 0; sx < 2; ++sx) {
                for (int sy = 0; sy < 2; ++sy) {
                    if (!owns_edge(vbx, dom, sx ? REMORASeaIce::IceSide::east
                                                : REMORASeaIce::IceSide::west)) { continue; }
                    if (!owns_edge(vbx, dom, sy ? REMORASeaIce::IceSide::north
                                                : REMORASeaIce::IceSide::south)) { continue; }
                    const int ci = cx[sx], cj = cy[sy];
                    const int di = dx[sx], dj = dx[sy];
                    Box corner(IntVect(ci,cj,0), IntVect(ci,cj,0));
                    ParallelFor(corner, [=] AMREX_GPU_DEVICE (int i, int j, int)
                    {
                        Real e, tt;
                        REMORASeaIce::ice_tibc_corner(
                            i, j, di, dj, Array4<Real const>(s_new),
                            Array4<Real const>(s_new), Array4<Real const>(s_new),
                            mskr, min_hi, e, tt,
                            REMORASeaIce::IceS::Hice, REMORASeaIce::IceS::Tice,
                            REMORASeaIce::IceS::Enth);
                        s_new(i,j,0,REMORASeaIce::IceS::Enth) = e;
                        s_new(i,j,0,REMORASeaIce::IceS::Tice) = tt;
                    });
                }
            }
        }
    }
}

/**
 * \brief The ice velocity boundary conditions: ROMS ice_uibc.F and ice_vibc.F.
 *
 * u is NORMAL on the west and east, where the boundary FACE value is set and
 * there is no ghost cell to fill, and TANGENTIAL on the south and north,
 * where a ghost row is filled like a scalar's. v is the other way round. The
 * two are not symmetric in ROMS either -- three of the four tangential
 * `closed` blocks scale the interior value and the fourth scales the ghost
 * cell in place -- so both components go through the same kernels with the
 * side telling them which case they are in.
 */
void
REMORA::ice_velbc (int lev, MultiFab& mf_u_new, MultiFab& mf_v_new,
                   const MultiFab& mf_u_old, const MultiFab& mf_v_old,
                   Real dtice)
{
    BL_PROFILE("REMORA::ice_velbc()");
    if (!REMORASeaIce::has_sea_ice(sea_ice_model)) { return; }

    const auto& p   = sea_ice_params;
    const Box& dom  = geom[lev].Domain();
    const auto& per = geom[lev].periodicity();
    const Real gamma2 = p.bc_gamma2;
    const auto fix = p.bc_fix;

    namespace IceS = REMORASeaIce::IceS;
    const int bc_u = REMORASeaIce::SeaIceParameters::ice_bc_var_for_comp(IceS::Uice);
    const int bc_v = REMORASeaIce::SeaIceParameters::ice_bc_var_for_comp(IceS::Vice);

    // Each component in turn, and within it ROMS's side order.
    for (int uv = 0; uv < 2; ++uv)
    {
        const bool do_u  = (uv == 0);
        const int  bc_var = do_u ? bc_u : bc_v;
        MultiFab& V_new       = do_u ? mf_u_new : mf_v_new;
        const MultiFab& V_old = do_u ? mf_u_old : mf_v_old;

        for (int side = 0; side < REMORASeaIce::IceSide::count; ++side)
        {
            const auto type = p.bc_type[bc_var][side];
            if (type == REMORASeaIce::IceBCType::periodic) { continue; }
            const bool side_is_x = REMORASeaIce::ice_side_is_x(side);
            if (per.isPeriodic(side_is_x ? 0 : 1)) { continue; }

            // Is this component NORMAL to this side? u is normal on an
            // x-normal side, v on a y-normal one.
            const bool normal = (do_u == side_is_x);
            const Real bval   = p.bc_value[bc_var][side];
            const bool nudge  = p.bc_nudging[side];
            const Real tau_in = p.bc_tau_in, tau_out = p.bc_tau_out;
            const bool closed = (type == REMORASeaIce::IceBCType::closed);
            const auto bdat = ice_bdry_array(lev, bc_var, side);
            const bool from_file = (bdat.p != nullptr);

            for (MFIter mfi(V_new, false); mfi.isValid(); ++mfi)
            {
                const Box vbx = convert(mfi.validbox(), IntVect::TheZeroVector());
                if (!owns_edge(vbx, dom, side)) { continue; }

                Array4<Real      > const& vn = V_new.array(mfi);
                Array4<Real const> const& vo = V_old.const_array(mfi);
                Array4<Real const> const& msk =
                    do_u ? vec_msku[lev]->const_array(mfi)
                         : vec_mskv[lev]->const_array(mfi);

                // The boundary point. For a velocity NORMAL to the side it is
                // the face itself -- Istr on the west, Iend+1 on the east --
                // and for a tangential one it is the last interior row, with
                // the ghost one step out.
                int ib = 0, jb = 0;
                switch (side) {
                    case REMORASeaIce::IceSide::west:
                        ib = dom.smallEnd(0); break;
                    case REMORASeaIce::IceSide::east:
                        ib = normal ? (dom.bigEnd(0) + 1) : dom.bigEnd(0); break;
                    case REMORASeaIce::IceSide::south:
                        jb = dom.smallEnd(1); break;
                    default:
                        jb = normal ? (dom.bigEnd(1) + 1) : dom.bigEnd(1); break;
                }
                // The sweep along the side. A tangential component has one
                // more point along the side than the box has cells. Where the
                // box meets another box or a periodic edge that point takes
                // the condition too: ROMS exchanges the boundary row after
                // the condition (ice_elastic.F:738-754), and FillBoundary
                // does not, the row being outside every box. At a physical
                // edge the point is the corner, which only ROMS's tangential
                // `closed` reaches, running over the TOTAL range so that free
                // slip fills the corner face as well.
                const int  t     = side_is_x ? 1 : 0;
                const bool t_per = per.isPeriodic(t);
                const int  tlo   = dom.smallEnd(t), tlen = dom.length(t);
                const int lo = std::max(vbx.smallEnd(t), dom.smallEnd(t));
                int hi = std::min(vbx.bigEnd(t), dom.bigEnd(t));
                if (!normal && (closed || t_per || vbx.bigEnd(t) < dom.bigEnd(t))) { hi += 1; }
                if (hi < lo) { continue; }
                Box line(IntVect(side_is_x ? ib : lo, side_is_x ? lo : jb, 0),
                         IntVect(side_is_x ? ib : hi, side_is_x ? hi : jb, 0));

                const int l_side = side;
                const auto l_type = type;
                if (normal) {
                    ParallelFor(line, [=] AMREX_GPU_DEVICE (int i, int j, int)
                    {
                        // The boundary FACE is the cell being written, and
                        // the fab's box sits on it for a u-type record on an
                        // x side (or a v-type on a y side), so no offset.
                        const Real bv = from_file ? bdat(i,j,0) : bval;
                        vn(i,j,0) = REMORASeaIce::ice_vel_bc_cell(
                            i, j, l_side, l_type, vo, Array4<Real const>(vn),
                            msk, nudge, bv, tau_in, tau_out, dtice, fix);
                    });
                } else {
                    ParallelFor(line, [=] AMREX_GPU_DEVICE (int i, int j, int)
                    {
                        const int d  = REMORASeaIce::ice_side_dir(l_side);
                        const int ig = side_is_x ? (i - d) : i;
                        const int jg = side_is_x ? j       : (j - d);
                        // The last point on a periodic side takes the
                        // value of the first, which it images.
                        Real bv = bval;
                        if (from_file) {
                            int tb = side_is_x ? jg : ig;
                            if (t_per) { tb = tlo + ((tb - tlo) % tlen + tlen) % tlen; }
                            bv = side_is_x ? bdat(ig,tb,0) : bdat(tb,jg,0);
                        }
                        vn(ig,jg,0) = REMORASeaIce::ice_vel_bc_tangential(
                            i, j, l_side, l_type, vo, Array4<Real const>(vn),
                            msk, nudge, bv, tau_in, tau_out, dtice,
                            gamma2, fix);
                    });
                }
            }
        }

        // The corners. On a u point the x corner index is the FACE (Istr or
        // Iend+1) and the y one the ghost row; on a v point it is the other
        // way round.
        if (per.isPeriodic(0) || per.isPeriodic(1)) { continue; }
        for (MFIter mfi(V_new, false); mfi.isValid(); ++mfi)
        {
            const Box vbx = convert(mfi.validbox(), IntVect::TheZeroVector());
            Array4<Real> const& vn = V_new.array(mfi);
            const int cx[2] = {do_u ? dom.smallEnd(0) : dom.smallEnd(0) - 1,
                               dom.bigEnd(0) + 1};
            const int cy[2] = {do_u ? dom.smallEnd(1) - 1 : dom.smallEnd(1),
                               dom.bigEnd(1) + 1};
            const int dx[2] = {+1, -1};
            for (int sx = 0; sx < 2; ++sx) {
                for (int sy = 0; sy < 2; ++sy) {
                    if (!owns_edge(vbx, dom, sx ? REMORASeaIce::IceSide::east
                                                : REMORASeaIce::IceSide::west)) { continue; }
                    if (!owns_edge(vbx, dom, sy ? REMORASeaIce::IceSide::north
                                                : REMORASeaIce::IceSide::south)) { continue; }
                    const int ci = cx[sx], cj = cy[sy];
                    const int di = dx[sx], dj = dx[sy];
                    Box corner(IntVect(ci,cj,0), IntVect(ci,cj,0));
                    ParallelFor(corner, [=] AMREX_GPU_DEVICE (int i, int j, int)
                    {
                        vn(i,j,0) = REMORASeaIce::ice_bc_corner(
                            i, j, di, dj, Array4<Real const>(vn));
                    });
                }
            }
        }
    }
}
