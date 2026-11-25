#ifndef REMORA_ICEBC_H_
#define REMORA_ICEBC_H_

#include <algorithm>
#include <cmath>

#include <AMReX_REAL.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>

#include <REMORA_SeaIce.H>

namespace REMORASeaIce {

/**
 * Lateral boundary conditions for the ice, ported from ROMS ice_bc2d.F,
 * ice_tibc.F, ice_uibc.F and ice_vibc.F.
 *
 * ROMS writes each of the four sides out separately in each of the four routines. Except for
 * the cases the kernels note, the sides differ only in
 *
 *   - which index is the boundary-adjacent interior cell (Istr/Iend/Jstr/Jend)
 *   - the direction to the second interior cell: +1 on a lo side, -1 on a hi
 *   - which ghost cell is written: ib-d
 *   - whether the normal is x or y
 *
 * so one implementation parameterised by (normal, d) covers all four. The old level is ROMS's
 * `know`. ROMS's RADIATION_2D is not ported, so the phase speed along the side is zero and the
 * terms it multiplies are left out.
 */
namespace IceBC {
    constexpr amrex::Real eps = amrex::Real(1.0e-20);
}

/** Is this side normal to x? */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
constexpr bool ice_side_is_x (const int side)
{
    return (side == IceSide::west) || (side == IceSide::east);
}

/** +1 on a lo side (west, south), -1 on a hi side (east, north). */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
constexpr int ice_side_dir (const int side)
{
    return ((side == IceSide::west) || (side == IceSide::south)) ? 1 : -1;
}

/**
 * The tangential difference ROMS calls `grad`, masked.
 *
 * On an x-normal side it differences in j and masks with vmask; on a
 * y-normal side it differences in i and masks with umask. `at` is the cell
 * the difference is centred on.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
ice_bc_grad (const int i, const int j, const bool is_x,
             const amrex::Array4<amrex::Real const>& f_know,
             const amrex::Array4<amrex::Real const>& umask,
             const amrex::Array4<amrex::Real const>& vmask,
             const int nc = 0)
{
    if (is_x) {
        return (f_know(i,j,0,nc) - f_know(i,j-1,0,nc)) * vmask(i,j,0);
    }
    return (f_know(i,j,0,nc) - f_know(i-1,j,0,nc)) * umask(i,j,0);
}

/**
 * One ghost cell of the radiation condition for an ice SCALAR at rho points.
 *
 * `ib`/`jb` is the boundary-adjacent INTERIOR cell. The ghost cell written is
 * one step outward from it.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
ice_bc_radiation (const int ib, const int jb, const int side,
                  const amrex::Array4<amrex::Real const>& f_know,
                  const amrex::Array4<amrex::Real const>& f_new,
                  const amrex::Array4<amrex::Real const>& umask,
                  const amrex::Array4<amrex::Real const>& vmask,
                  const amrex::Array4<amrex::Real const>& rmask,
                  const bool nudging,
                  const amrex::Real bval,       //!< boundary value for nudging
                  const amrex::Real tobc_in,
                  const amrex::Real tobc_out,
                  const amrex::Real dt,
                  //! Which component of a multi-component ice state this is.
                  //! The masks are always single-component.
                  const int nc = 0)
{
    using amrex::Real;
    const bool is_x = ice_side_is_x(side);
    const int  d    = ice_side_dir(side);

    // The ghost cell and the second interior cell.
    const int ig = is_x ? (ib - d) : ib;
    const int jg = is_x ? jb       : (jb - d);
    const int i2 = is_x ? (ib + d) : ib;
    const int j2 = is_x ? jb       : (jb + d);
    // The tangential neighbour, +1 along the boundary.
    const int it = is_x ? 0 : 1;   // 1 => step in i, 0 => step in j
    const int ib_t = ib + it;
    const int jb_t = jb + (1 - it);

    const Real dTdt = f_know(ib,jb,0,nc) - f_new(ib,jb,0,nc);
    // The NORMAL gradient: dTdx on an x side, dTde on a y side. ROMS names
    // them differently per side but the expression is the same.
    const Real dTdn = f_new(ib,jb,0,nc) - f_new(i2,j2,0,nc);

    Real tau = Real(0.0);
    if (nudging) {
        tau = tobc_out;
        if ((dTdt * dTdn) < Real(0.0)) { tau = tobc_in; }
        tau *= dt;
    }
    const Real dTdt_lim = ((dTdt * dTdn) < Real(0.0)) ? Real(0.0) : dTdt;

    // The TANGENTIAL gradient, upstream-selected on the sign of dTdt.
    const Real g_b   = ice_bc_grad(ib,   jb,   is_x, f_know, umask, vmask, nc);
    const Real g_b_t = ice_bc_grad(ib_t, jb_t, is_x, f_know, umask, vmask, nc);
    const Real dTdt_ = dTdt_lim;
    const Real dTdtg = dTdt_ * (g_b + g_b_t);
    const Real dTdtan = (dTdtg > Real(0.0)) ? g_b : g_b_t;

    const Real cff = amrex::max(dTdn * dTdn + dTdtan * dTdtan, Real(IceBC::eps));
    const Real Cn = dTdt_ * dTdn;

    Real out = (cff * f_know(ig,jg,0,nc) + Cn * f_new(ib,jb,0,nc)) / (cff + Cn);
    if (nudging) {
        out += tau * (bval - f_know(ig,jg,0,nc));
    }
    return out * rmask(ig,jg,0);
}

/**
 * The non-radiation conditions for an ice scalar: clamped, mixed and closed.
 *
 * `vel_inward` is the ice velocity normal to the side at the boundary FACE,
 * signed so that positive means into the domain. `mixed` uses it to pick
 * inflow from outflow.
 *
 * ROMS is inconsistent about that index and it is worth knowing: the lo sides
 * test a hard-coded `ui(1,j)` and `vi(i,1)` while the hi sides test the
 * symbolic `ui(Iend+1,j)` and `vi(i,Jend+1)` (ice_bc2d.F:150, 380 against
 * :265, 495). The literal is correct only because `DOMAIN%Western_Edge(tile)`
 * restricts the block to the tile that owns index 1, where Istr == 1. The
 * port uses the boundary face symbolically on all four sides, which is what
 * ROMS means; on a single-tile global grid the two agree exactly.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
ice_bc_simple (const int ib, const int jb, const int side,
               const IceBCType type,
               const amrex::Array4<amrex::Real const>& f_old,
               const amrex::Array4<amrex::Real const>& f_new,
               const amrex::Array4<amrex::Real const>& rmask,
               const amrex::Real bval,
               const amrex::Real vel_inward,
               const int nc = 0)
{
    using amrex::Real;
    const bool is_x = ice_side_is_x(side);
    const int  d    = ice_side_dir(side);
    const int ig = is_x ? (ib - d) : ib;
    const int jg = is_x ? jb       : (jb - d);

    Real out;
    switch (type) {
        case IceBCType::clamped:
            out = bval;
            break;
        case IceBCType::mixed:
            // Inflow takes the boundary value; outflow takes the OLD interior
            // value, which is upstream of the face.
            out = (vel_inward >= Real(0.0)) ? bval : f_old(ib,jb,0,nc);
            break;
        case IceBCType::closed:
        case IceBCType::gradient:
            // ROMS's ice `closed` IS a gradient condition, not a wall.
            out = f_new(ib,jb,0,nc);
            break;
        default:
            return f_new(ig,jg,0,nc);   // none / periodic: leave it alone
    }
    return out * rmask(ig,jg,0);
}

/**
 * The ice-temperature boundary condition, ported from ROMS ice_tibc.F.
 *
 * This is NOT ice_bc2d applied to Tice. ROMS carries the boundary condition on
 * the ENTHALPY, `hi*ti`, and recovers the temperature from it -- so a
 * clamped side takes the PRODUCT of the boundary thickness and temperature,
 * and both the enthalpy and the temperature are zeroed wherever the ghost
 * thickness is at or below min_hi. There is no radiation branch.
 *
 * `gradient` and `closed` are NOT the same here, and look it only at a
 * glance: both set the ghost enthalpy from the interior, but `gradient` then
 * derives `ti` from the GHOST thickness while `closed` copies `ti` straight
 * from the interior. They differ whenever hi(ghost) != hi(interior).
 *
 * ROMS's `mixed` is not the same on all four sides. The temperature is recovered from the ghost
 * enthalpy and the min_hi guard applied
 *
 *   - OUTSIDE the inflow test on the west and the north
 *    , so both branches get a temperature and the guard
 *   - INSIDE the outflow branch on the east and the south, so
 *     an east or south side with INFLOW sets the ghost enthalpy from the
 *     boundary value, leaves `ti` at whatever the ghost cell already held,
 *     and never applies the guard at all
 *
 * The four are otherwise the same condition, and three of them write the
 * outflow enthalpy as the stored `enthalpy(interior,liold)` while the east
 * writes the product `hi(interior,liold)*ti(interior,liold)` -- equal only
 * because ice_enthalpy.F keeps the two consistent.
 *
 * With `remora.ice.bc_fix.tibc_mixed`, all four sides behave like the west and the north.
 *
 * Two more places where one side of the four is written differently, both
 * with no effect while the enthalpy array agrees with hi*ti -- which
 * ice_enthalpy.F maintains -- and both reproduced anyway so that the
 * comparison against ROMS holds even when they do not agree:
 *
 *   - `closed` takes the ghost enthalpy from hi*ti at the new level on the
 *     west, east and north, but from the stored `enthalpy` array on the
 *     south
 *   - `mixed` on OUTFLOW takes it from the stored `enthalpy` at the old level
 *     on the west, south and north, but from the product hi*ti at the old
 *     level on the east
 *
 * So the kernel reads all six arrays ROMS reads: thickness, temperature and
 * enthalpy at each of the two levels. Nothing is inferred from another.
 *
 * One CPP branch is NOT taken and is recorded so it is not mistaken for an
 * omission: under `SOGLOBEC` the northern `gradient` divides by the INTERIOR
 * thickness rather than the ghost one. The port follows the
 * #else, like the other sides.
 *
 * Writes both the ghost enthalpy and the ghost temperature. When the stale-ti
 * case above applies it writes the enthalpy ONLY, so the caller must pass
 * `ti_ghost` already holding the ghost cell's current temperature -- which is
 * what makes "leave it alone" expressible at all.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
ice_tibc_cell (const int ib, const int jb, const int side,
               const IceBCType type,
               const amrex::Array4<amrex::Real const>& hi_old,
               const amrex::Array4<amrex::Real const>& ti_old,
               const amrex::Array4<amrex::Real const>& hi_new,
               const amrex::Array4<amrex::Real const>& ti_new,
               const amrex::Array4<amrex::Real const>& enth_old,
               const amrex::Array4<amrex::Real const>& enth_new,
               const amrex::Array4<amrex::Real const>& rmask,
               const amrex::Real bval_hice,   //!< boundary ice thickness
               const amrex::Real bval_tice,   //!< boundary ice temperature
               const amrex::Real vel_inward,
               const amrex::Real min_hi,
               const IceBCQuirkFix& fix,
               amrex::Real& enth_ghost,
               amrex::Real& ti_ghost,
               const int c_hi = 0, const int c_ti = 0, const int c_enth = 0)
{
    using amrex::Real;
    const bool is_x = ice_side_is_x(side);
    const int  d    = ice_side_dir(side);
    const int ig = is_x ? (ib - d) : ib;
    const int jg = is_x ? jb       : (jb - d);
    const Real m  = rmask(ig,jg,0);
    // ROMS's eps in ice_tibc.F, the floor on the thickness in that division.
    constexpr Real eps_t = Real(1.0e-6);

    bool ti_from_interior = false;
    switch (type) {
        case IceBCType::clamped:
            enth_ghost = bval_hice * bval_tice * m;
            break;
        case IceBCType::mixed: {
            const bool inflow = (vel_inward >= Real(0.0));
            // The east writes the product where the other three write the
            // stored array; both are the old level.
            const Real outflow_enth = (side == IceSide::east)
                ? (hi_old(ib,jb,0,c_hi) * ti_old(ib,jb,0,c_ti))
                : enth_old(ib,jb,0,c_enth);
            enth_ghost = (inflow ? (bval_hice * bval_tice)
                                 : outflow_enth) * m;
            // East and south put the ti recovery and the min_hi guard inside
            // the outflow branch, so an INFLOW there leaves both alone.
            const bool stale = inflow && !fix.tibc_mixed
                            && ((side == IceSide::east) ||
                                (side == IceSide::south));
            if (stale) { return; }
            break;
        }
        case IceBCType::gradient:
            enth_ghost = hi_new(ib,jb,0,c_hi) * ti_new(ib,jb,0,c_ti) * m;
            break;
        case IceBCType::closed:
            // The south reads the stored array where the other three form
            // the product; both are the new level.
            enth_ghost = ((side == IceSide::south)
                          ? enth_new(ib,jb,0,c_enth)
                          : (hi_new(ib,jb,0,c_hi) * ti_new(ib,jb,0,c_ti))) * m;
            ti_from_interior = true;
            break;
        default:
            return;   // none / periodic
    }

    ti_ghost = ti_from_interior
        ? (ti_new(ib,jb,0,c_ti) * m)
        : (enth_ghost / amrex::max(hi_new(ig,jg,0,c_hi), eps_t));

    if (hi_new(ig,jg,0,c_hi) <= min_hi) {
        enth_ghost = Real(0.0);
        ti_ghost   = Real(0.0);
    }
}

/**
 * The ice-velocity boundary condition, ported from ROMS ice_uibc.F and
 * ice_vibc.F.
 *
 * The difference from the scalar case is the staggering: the velocity point
 * sits ON the boundary face, so there is no ghost cell to fill -- the
 * boundary VALUE itself is set, and the radiation reads one cell further in.
 * For the western edge ROMS writes `ui(Istr,j)` using `ui(Istr+1,j)` as the
 * interior and `ui(Istr+2,j)` beyond it, where the scalar wrote the ghost
 * `field(Istr-1,j)` from `field(Istr,j)`.
 *
 * `ib`/`jb` is the boundary velocity point; one and two points inward are
 * derived from (side, d) as everywhere else. `mask` is the mask AT that
 * velocity point -- umask for u, vmask for v.
 *
 * Two details that differ from the scalar case:
 *
 *   - the tangential `grad` is NOT masked here, where ice_bc2d.F multiplies
 *     its own by vmask
 *   - `closed` sets the normal velocity to exactly ZERO -- a wall -- and
 *     applies NO mask, where `gradient` copies the neighbour and does mask
 *
 * Every one of the eight normal blocks across the two files masks with the
 * mask AT the point being written, except the northern `clamped` and `gradient` of
 * ice_vibc.F, which use the one a cell inside. `remora.ice.bc_fix.vibc_north_mask` uses the
 * boundary point there too.
 *
 * ice_vibc.F's southern radiation block is guarded by the condition of another variable's
 * record (`LBC(isouth,isVice,ng)` for `ibICE(isVice)`); here it is Vice's own.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
ice_vel_bc_cell (const int ib, const int jb, const int side,
                 const IceBCType type,
                 const amrex::Array4<amrex::Real const>& v_know,
                 const amrex::Array4<amrex::Real const>& v_new,
                 const amrex::Array4<amrex::Real const>& mask,
                 const bool nudging,
                 const amrex::Real bval,
                 const amrex::Real tau_in,
                 const amrex::Real tau_out,
                 const amrex::Real dt,
                 const IceBCQuirkFix& fix)
{
    using amrex::Real;
    const bool is_x = ice_side_is_x(side);
    const int  d    = ice_side_dir(side);
    // One and two points inward from the boundary face.
    const int i1 = is_x ? (ib + d) : ib;
    const int j1 = is_x ? jb       : (jb + d);
    // The mask these two ROMS blocks reach for; see the note above. `north`
    // on a normal side can only be the v component.
    const bool mask_inward = (side == IceSide::north) && !fix.vibc_north_mask;
    const int  im = mask_inward ? i1 : ib;
    const int  jm = mask_inward ? j1 : jb;
    const int i2 = is_x ? (ib + 2*d) : ib;
    const int j2 = is_x ? jb         : (jb + 2*d);
    const int it   = is_x ? 0 : 1;
    const int i1_t = i1 + it, j1_t = j1 + (1 - it);

    switch (type) {
        case IceBCType::clamped:
            return bval * mask(im,jm,0);
        case IceBCType::gradient:
        case IceBCType::closed:
            // ROMS's closed on a velocity point is a ZERO normal flow, not a
            // copy: a wall. gradient copies the neighbour.
            return (type == IceBCType::closed)
                 ? Real(0.0)
                 : v_new(i1,j1,0) * mask(im,jm,0);
        case IceBCType::radiation: break;
        default: return v_new(ib,jb,0);
    }

    const Real dUdt = v_know(i1,j1,0) - v_new(i1,j1,0);
    const Real dUdn = v_new(i1,j1,0)  - v_new(i2,j2,0);

    Real tau = Real(0.0);
    if (nudging) {
        tau = ((dUdt * dUdn) < Real(0.0)) ? tau_in : tau_out;
        tau *= dt;
    }
    const Real dUdt_lim = ((dUdt * dUdn) < Real(0.0)) ? Real(0.0) : dUdt;

    const Real g1   = v_know(i1,  j1,  0) - v_know(i1  - it, j1 - (1-it), 0);
    const Real g1_t = v_know(i1_t,j1_t,0) - v_know(i1_t - it, j1_t - (1-it), 0);
    const Real dUtan = ((dUdt_lim * (g1 + g1_t)) > Real(0.0)) ? g1 : g1_t;

    const Real cff = amrex::max(dUdn * dUdn + dUtan * dUtan, Real(IceBC::eps));
    const Real Cn  = dUdt_lim * dUdn;

    Real out = (cff * v_know(ib,jb,0) + Cn * v_new(i1,j1,0)) / (cff + Cn);
    if (nudging) {
        out += tau * (bval - v_know(ib,jb,0));
    }
    return out * mask(ib,jb,0);
}

/**
 * The ice velocity on a side where it is TANGENTIAL, from ROMS ice_uibc.F's
 * southern/northern blocks and ice_vibc.F's western/eastern ones.
 *
 * This is a different condition from the normal case and behaves like a
 * scalar: it fills a GHOST row rather than setting a boundary point, because
 * a u point is interior to a y-normal boundary. Three things differ from both
 * the normal-velocity and the scalar cases:
 *
 *   - the tangential `grad` is a FORWARD difference, `u(i+1) - u(i)`, where
 *     the scalar's is backward, so the upstream pick brackets the point with
 *     `grad(i-1)` and `grad(i)` instead of `grad(i)` and `grad(i+1)`
 *   - `closed` is `gamma2 * u(interior)` -- free slip at gamma2 = +1, no slip
 *     at -1 -- not zero and not a copy
 *   - `gradient` copies the interior, and both are masked
 *
 * ROMS writes the interior index as a literal `1` in some `gradient` and `closed` branches,
 * correct only on the tile that owns index 1. The port uses the symbolic index.
 *
 * On the northern side of ice_uibc.F the nudging is inside a CPP guard that nothing defines,
 * so that side does not nudge. `remora.ice.bc_fix.uibc_north_nudging` makes it. A call with
 * `side == north` can only be the u component.
 *
 * Three of the four tangential `closed` blocks scale the adjacent interior value by gamma2,
 * which is free slip at +1 and no slip at -1. The western v block scales the ghost cell in
 * place. `remora.ice.bc_fix.vibc_west_closed` gives it the form of the others. A call with
 * `side == west` can only be the v component.
 *
 * `ib`/`jb` is the boundary-adjacent INTERIOR point; the ghost filled is one
 * step outward. `tan_lo`/`tan_hi` are that point's tangential neighbours, so
 * the caller does not have to know which axis is tangential.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
ice_vel_bc_tangential (const int ib, const int jb, const int side,
                       const IceBCType type,
                       const amrex::Array4<amrex::Real const>& v_know,
                       const amrex::Array4<amrex::Real const>& v_new,
                       const amrex::Array4<amrex::Real const>& mask,
                       const bool nudging,
                       const amrex::Real bval,
                       const amrex::Real tau_in,
                       const amrex::Real tau_out,
                       const amrex::Real dt,
                       const amrex::Real gamma2,
                       const IceBCQuirkFix& fix)
{
    using amrex::Real;
    // On a y-normal side the u component is tangential, so the NORMAL
    // direction here is y and the tangential one is x -- the opposite of
    // ice_side_is_x's meaning for the component being set.
    const bool norm_is_x = ice_side_is_x(side);
    const int  d = ice_side_dir(side);
    const int ig = norm_is_x ? (ib - d) : ib;
    const int jg = norm_is_x ? jb       : (jb - d);
    const int i2 = norm_is_x ? (ib + d) : ib;
    const int j2 = norm_is_x ? jb       : (jb + d);
    // One step along the boundary.
    const int ts = norm_is_x ? 0 : 1;

    switch (type) {
        case IceBCType::clamped:
            return bval * mask(ig,jg,0);
        case IceBCType::gradient:
            return v_new(ib,jb,0) * mask(ig,jg,0);
        case IceBCType::closed: {
            // Free slip at +1, no slip at -1. NOT zero: the normal component
            // is the one that vanishes at a wall, not the tangential.
            // ice_vibc.F's WEST scales the ghost cell in place rather than
            // the interior -- see the note above.
            const bool in_place = (side == IceSide::west)
                               && !fix.vibc_west_closed;
            return gamma2 * (in_place ? v_new(ig,jg,0) : v_new(ib,jb,0))
                 * mask(ig,jg,0);
        }
        case IceBCType::radiation: break;
        default: return v_new(ig,jg,0);
    }

    const Real dUdt = v_know(ib,jb,0) - v_new(ib,jb,0);
    const Real dUdn = v_new(ib,jb,0)  - v_new(i2,j2,0);

    // See the note above on the northern side
    const bool nudge = nudging
                    && ((side != IceSide::north) || fix.uibc_north_nudging);
    Real tau = Real(0.0);
    if (nudge) {
        tau = ((dUdt * dUdn) < Real(0.0)) ? tau_in : tau_out;
        tau *= dt;
    }
    const Real dUdt_lim = ((dUdt * dUdn) < Real(0.0)) ? Real(0.0) : dUdt;

    // Forward differences along the boundary, at the interior row and at the
    // ghost row, with the pair that brackets the point being set.
    auto fwd = [&] (int i, int j) {
        return v_know(i + ts, j + (1 - ts), 0) - v_know(i, j, 0);
    };
    const Real gb_lo = fwd(ib - ts, jb - (1 - ts));
    const Real gb_hi = fwd(ib, jb);
    const Real dUtan = ((dUdt_lim * (gb_lo + gb_hi)) > Real(0.0)) ? gb_lo : gb_hi;

    const Real cff = amrex::max(dUdn * dUdn + dUtan * dUtan, Real(IceBC::eps));
    const Real Cn  = dUdt_lim * dUdn;

    Real out = (cff * v_know(ig,jg,0) + Cn * v_new(ib,jb,0)) / (cff + Cn);
    if (nudge) {
        out += tau * (bval - v_know(ig,jg,0));
    }
    return out * mask(ig,jg,0);
}

/**
 * Boundary corners, from the block at the tail of ROMS's ice_bc2d.F,
 * ice_tibc.F, ice_uibc.F and ice_vibc.F.
 *
 * The corner cell is the mean of its two neighbours one step inward, one along x and one
 * along y.
 *
 * The caller supplies the corner cell and the two inward steps, because WHICH
 * cell is the corner depends on the stagger -- on the western side a rho or
 * v point has its ghost column at Istr-1 while a u point has its boundary
 * face at Istr, so for a south-west corner:
 *
 *   | field | corner cell       | di | dj |
 *   |-------|-------------------|----|----|
 *   | rho   | (Istr-1, Jstr-1)  | +1 | +1 |
 *   | u     | (Istr  , Jstr-1)  | +1 | +1 |
 *   | v     | (Istr-1, Jstr  )  | +1 | +1 |
 *
 * and di/dj flip to -1 on the eastern/northern side, with the corner cell at
 * Iend+1/Jend+1 for every stagger there.
 *
 * ROMS guards the whole block with `.not.(EWperiodic .or. NSperiodic)`, so a
 * periodic direction suppresses all four corners, not just the two on that
 * axis. ice_bc2d.F spells the same guard `.not.EWperiodic .and. .not.
 * NSperiodic`, which is the same thing by De Morgan.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
ice_bc_corner (const int i, const int j, const int di, const int dj,
               const amrex::Array4<amrex::Real const>& f,
               const int nc = 0)
{
    return amrex::Real(0.5) * (f(i+di,j,0,nc) + f(i,j+dj,0,nc));
}

/**
 * The same corner for the ice temperature, from ice_tibc.F's corner block.
 *
 * Two things it does that the plain corner does not, both ROMS's:
 *
 *   - `ti` is AVERAGED from its own two neighbours, not divided out of the
 *     averaged enthalpy. ROMS carries the division as a commented-out line
 *     immediately above, so the averaging
 *     replaced it deliberately and the port follows what runs.
 *   - both fields are masked with rmask at the corner and then zeroed
 *     outright when the corner thickness is at or below min_hi, which is the
 *     same guard ice_tibc_cell applies on the sides.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
ice_tibc_corner (const int i, const int j, const int di, const int dj,
                 const amrex::Array4<amrex::Real const>& enth_new,
                 const amrex::Array4<amrex::Real const>& ti_new,
                 const amrex::Array4<amrex::Real const>& hi_new,
                 const amrex::Array4<amrex::Real const>& rmask,
                 const amrex::Real min_hi,
                 amrex::Real& enth_corner,
                 amrex::Real& ti_corner,
                 const int c_hi = 0, const int c_ti = 0, const int c_enth = 0)
{
    using amrex::Real;
    const Real m = rmask(i,j,0);
    enth_corner = ice_bc_corner(i, j, di, dj, enth_new, c_enth) * m;
    ti_corner   = ice_bc_corner(i, j, di, dj, ti_new,   c_ti)   * m;
    if (hi_new(i,j,0,c_hi) <= min_hi) {
        enth_corner = Real(0.0);
        ti_corner   = Real(0.0);
    }
}

} // namespace REMORASeaIce

#endif
