#ifndef REMORA_SEAICE_H_
#define REMORA_SEAICE_H_

#include <string>
#include <type_traits>

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

/**
 * Sea-ice model, ported from ROMS ICE_MODEL (ROMS/Nonlinear/SeaIce/).
 *
 * Budgell (2005) elastic-viscous-plastic dynamics with Mellor-Kantha (1989)
 * thermodynamics. Docs/sphinx_doc/SeaIce.rst describes the model, its options
 * and every departure from ROMS.
 */
namespace REMORASeaIce {

/**
 * Which sea-ice physics to run.
 *
 * One formulation, as in ROMS: Mellor and Kantha thermodynamics with
 * elastic-viscous-plastic momentum. Which parts of it run is set by the
 * `remora.ice.do_thermo`, `remora.ice.do_momentum` and `remora.ice.do_advect` switches, which
 * correspond to ROMS's ICE_THERMO, ICE_MOMENTUM and ICE_ADVECT.
 */
enum class SeaIceModel {
    none,          //!< default; reproduces a run without ice exactly
    mellor_kantha  //!< the ROMS sea-ice model
};

inline bool has_sea_ice (SeaIceModel m) { return m != SeaIceModel::none; }

SeaIceModel parse_sea_ice_model (const std::string& name);
std::string sea_ice_model_name (SeaIceModel m);

/**
 * Component indices into the prognostic ice state, matching ROMS `Si` in
 * ice_mod.h **exactly, including the numbering**.
 *
 * The numbering is ROMS's, so that the two can be compared field by field. The cost is that the MultiFab carries a few components that hold nothing at
 * rho points -- `Uice`/`Vice`/`Uevp`/`Vevp` are at u and v points, so those
 * live in their own staggered MultiFabs and their slots here are never
 * written.
 *
 * `ISxy` is at rho points, not psi points: ice_evp_sig.F assigns
 * Si(i,j,...,isISxy) in an Istr..Iend, Jstr..Jend loop, and ice_elastic.F
 * averages four of them to reach a velocity point.
 *
 * `ice_biology` (isIphy, isINO3, isINH4, isIlog in ROMS) is not ported.
 */
namespace IceS {
    enum {
        Aice = 0,  //!< ice concentration [fraction]
        Hice,      //!< ice thickness, cell average [m]
        Hmel,      //!< melt-water thickness on ice [m]
        Hsno,      //!< snow thickness [m]
        Iage,      //!< ice age [days]; ROMS labels it seconds but adds dtice*sec2day
        ISxx,      //!< internal ice xx-stress
        ISxy,      //!< internal ice xy-stress, rho points
        ISyy,      //!< internal ice yy-stress
        Tice,      //!< ice interior temperature [C]
        Uice,      //!< ice U-velocity -- u points, unused here
        Vice,      //!< ice V-velocity -- v points, unused here
        Enth,      //!< ice/brine enthalpy
        Hage,      //!< thickness linked with ice age [m]
        Uevp,      //!< EVP ice U-velocity -- u points, unused here
        Vevp,      //!< EVP ice V-velocity -- v points, unused here
        NumComps
    };
}

/**
 * Component indices into the internal/diagnostic ice fields, matching ROMS
 * `Fi` in ice_mod.h exactly. Same reasoning as IceS: `AIus`/`AIvs` and
 * `Uavg`/`Vavg` are staggered and their slots here are unused.
 */
namespace IceF {
    enum {
        AIus = 0,  //!< surface air-ice U-stress -- u points, unused here
        AIvs,      //!< surface air-ice V-stress -- v points, unused here
        Bvis,      //!< ice bulk viscosity
        Hsse,      //!< sea surface elevation
        IOfv,      //!< ice-ocean friction velocity
        IOmf,      //!< ice-ocean mass flux
        IOmt,      //!< ice-ocean momentum transfer
        IOvs,      //!< ice-ocean velocity shear
        Isst,      //!< ice/snow surface temperature [C]
        Pgrd,      //!< gridded ice strength
        Pice,      //!< ice pressure or strength
        Qcon,      //!< ice/snow heat conductivity
        Qrhs,      //!< RHS heat flux over ice/snow
        Svis,      //!< ice shear viscosity
        S0mk,      //!< molecular sublayer salinity [psu]
        T0mk,      //!< molecular sublayer temperature [C]
        Uavg,      //!< mixed-layer U-velocity -- u points, unused here
        Vavg,      //!< mixed-layer V-velocity -- v points, unused here
        Wdiv,      //!< thickness change by advection per unit time [m/s]; ROMS's "divergence rate"
        W_ai,      //!< melt/freeze rate at air/ice [m/s]
        W_ao,      //!< melt/freeze rate at air/ocean [m/s]
        W_fr,      //!< ice accretion rate by frazil [m/s]
        W_io,      //!< melt/freeze rate at ice/ocean [m/s]
        W_ro,      //!< melt/freeze rate runoff [m/s]
        NumComps
    };
}

/**
 * Four places where ROMS's ice boundary routines write one side differently from the
 * other three in a way that changes the answer. The default, false, reproduces ROMS; true
 * makes that side behave like the others. The kernels are in REMORA_IceBC.H, and the keys
 * are `remora.ice.bc_fix.<member>`, or `remora.ice.bc_fix_all` for all of them.
 */
struct IceBCQuirkFix {
    //! Recover ti and apply the min_hi guard on an east or south `mixed`
    //! INFLOW, as the west and north blocks already do.
    bool tibc_mixed = false;
    //! Nudge the u velocity on the northern side, which ROMS leaves behind a
    //! never-defined CPP guard that would not compile if defined.
    bool uibc_north_nudging = false;
    //! Scale the INTERIOR value in the western v `closed`, not the ghost cell
    //! in place, so that gamma2 means free or no slip there as elsewhere.
    bool vibc_west_closed = false;
    //! Mask the northern v `clamped` and `gradient` with the mask AT the
    //! boundary point, not the one a cell inside it.
    bool vibc_north_mask = false;
};

/**
 * Which lateral condition an ice boundary side uses, under ROMS's LBC names.
 *
 * ROMS's ice `closed` is a gradient condition, not a wall: ice_bc2d.F has no separate
 * `gradient` branch, while ice_tibc.F, ice_uibc.F and ice_vibc.F have both and treat them
 * differently. Nudging is not a type: it is a switch that modifies `radiation`.
 */
enum class IceBCType {
    closed = 0, //!< ROMS %closed
    gradient,   //!< ROMS %gradient (distinct from closed outside ice_bc2d.F)
    radiation,  //!< ROMS %radiation, Orlanski-style
    clamped,    //!< ROMS %clamped: take the boundary value
    mixed,      //!< ROMS %mixed: clamped on inflow, upstream on outflow
    periodic    //!< handled by FillBoundary, so a no-op here
};

/** The four sides, in ROMS's own order: west, south, east, north. */
namespace IceSide { enum { west = 0, south, east, north, count }; }

IceBCType parse_ice_bc_type (const std::string& name);
std::string ice_bc_type_name (IceBCType t);

/**
 * Runtime parameters, read from `remora.ice.*`.
 *
 * Defaults are ROMS's own, from ROMS/External/ice.in. Where ROMS's file
 * carries a caveat in a comment it is repeated here, because the value and
 * the caveat belong together.
 */
struct SeaIceParameters {
    //! Number of elastic subcycles per ice timestep (ROMS nEVP).
    int nEVP = 60;

    //! Air density [kg/m3] (ROMS AirRho).
    amrex::Real AirRho = amrex::Real(1.267);

    //! Sea-ice density [kg/m3] (ROMS IceRho).
    amrex::Real IceRho = amrex::Real(900.0);

    //! Dry snow density [kg/m3] (ROMS SnowDryRho).
    amrex::Real SnowDryRho = amrex::Real(330.0);

    //! Wet snow density [kg/m3] (ROMS SnowWetRho).
    amrex::Real SnowWetRho = amrex::Real(450.0);

    //! Air-ice bulk drag coefficient (ROMS Cd_ai).
    amrex::Real Cd_ai = amrex::Real(3.0e-3);

    //! Ice-ocean bulk drag coefficient (ROMS Cd_io).
    amrex::Real Cd_io = amrex::Real(7.5e-3);

    //! Ice-strength exponential weighting on concentration (ROMS Astrength).
    amrex::Real Astrength = amrex::Real(20.0);

    //! Ice strength scale [N/m2] (ROMS Pstar).
    amrex::Real Pstar = amrex::Real(25000.0);

    // Minimum concentration and thickness limiters. ROMS/External/ice.in
    // carries the comment "I ran this with 0.002...think about" against both,
    // so these zeros are ROMS's shipped default and not a settled choice.
    //! Minimum fractional ice concentration (ROMS min_ai).
    amrex::Real min_ai = amrex::Real(0.0);

    //! Minimum average ice thickness [m] (ROMS min_hi).
    amrex::Real min_hi = amrex::Real(0.0);

    //! Maximum fractional ice concentration (ROMS max_ai).
    amrex::Real max_ai = amrex::Real(1.0);

    //! Maximum surface melt-water thickness [m] (ROMS max_hmelt).
    amrex::Real max_hmelt = amrex::Real(0.10);

    //! Minimum bulk viscosity [N s/m] (ROMS zetaMin).
    amrex::Real zetaMin = amrex::Real(0.0);

    //! Maximum bulk viscosity [N s/m] (ROMS zetaMax).
    amrex::Real zetaMax = amrex::Real(1.0e14);

    //! Turning angle for ice-water drag [rad] (ROMS stressAng).
    amrex::Real stressAng = amrex::Real(0.0);

    //! Squared ellipticity of the yield curve (ROMS ellip_sq).
    amrex::Real ellip_sq = amrex::Real(4.0);

    /**
     * Apply ROMS's own bulk-viscosity ceiling zmax = 2.5e8*Pice.
     *
     * false (default) reproduces ROMS, which computes zmax and then never
     * reads it, clamping with zetaMax (1e14 by default) instead -- leaving
     * the viscosity 20x to 2000x stiffer than zmax allows wherever the ice is
     * near-rigid. See ice_evp_viscosity() in REMORA_IceEVP.H.
     */
    bool zmax_clamp = false;

    /**
     * Include the pressure term in the frazil freezing point.
     *
     * false (default) reproduces ROMS, whose freezing_point(S, Z) takes a
     * depth and ignores it -- the live line is FP = -0.0543*S with both the
     * Gill (1982) and Steele (1989) forms commented out. The omitted term is
     * -0.076 C at 100 m, which is not small for a supercooling process. See
     * freezing_point() in REMORA_IceFrazil.H.
     */
    bool pressure_freezing = false;

    //! The departures from ROMS in the boundary routines; see IceBCQuirkFix.
    IceBCQuirkFix bc_fix;

    /**
     * Run the frazil-ice calculation (ROMS ice_frazil.F).
     *
     * Off by default because ROMS never calls ice_frazil: its only call is
     * inside ice_frazil.F itself, so frazil growth is always zero in a ROMS
     * run. It is the one part of the ice model that modifies the 3D tracers.
     */
    bool do_frazil = false;

    /**
     * Form Qnet_ao, the open-water heat flux the thermodynamics weights by
     * (1 - Aice), from the unattenuated shortwave.
     *
     * ROMS bulk_flux.F attenuates srflx by the ice cover first and forms
     * Qnet_ao from the result, so the default (false) does the same. True
     * uses the shortwave open water actually receives.
     */
    bool qnet_ao_open_water_sw = false;

    /**
     * ROMS ICE_MOM_BULK: bulk ice-ocean drag. The transfer coefficient is
     * Cd_io times the ice-ocean shear, in place of the log-layer form that
     * depends on ice thickness, and the air-ice stress on the ice is not
     * scaled by the thickness-dependent form-drag factor. Off unless the
     * ROMS application defines it; LAKE_ICE does.
     */
    bool mom_bulk = false;

    /**
     * ROMS ICE_CONVSNOW: convert to ice the snow whose base has been pushed
     * below sea level by its own weight. Off unless the ROMS application
     * defines it; LAKE_ICE does.
     */
    bool convsnow = false;

    /**
     * ROMS ICE_UPWIND: donor-cell upwind advection of the ice tracers, with
     * no antidiffusive correction. Off unless the ROMS application defines
     * it; LAKE_ICE does.
     */
    bool advect_upwind = false;

    /**
     * ROMS ICE_STRENGTH_QUAD: ice strength from the gridded parameter
     * Fi(icPgrd) and the square of the thickness, in place of Pstar and the
     * thickness. ROMS never fills that parameter, so this gives ice with no
     * strength; it is here so a ROMS run built with the option can be matched.
     */
    bool strength_quad = false;

    /**
     * Run the Mellor-Kantha thermodynamics (ROMS ICE_THERMO / ICE_MK). With
     * it off the ice state is carried by the dynamics alone and the ocean
     * keeps the surface fluxes the bulk flux computed.
     */
    bool do_thermo = true;

    /**
     * Run the EVP rheology and the ice momentum equations (ROMS
     * ICE_MOMENTUM). With it off the ice velocity stays as it was
     * initialised, which is ROMS's locally thermodynamic ice.
     */
    bool do_momentum = true;

    /** Advect the ice (ROMS ICE_ADVECT). */
    bool do_advect = true;

    /**
     * What happens to the ice state when it is not advected.
     *
     * ROMS keeps two copies of the ice state and switches which one is
     * current at the top of every ice step. The advection is the only thing
     * that fills the new current copy from the old one. Built without
     * ICE_ADVECT, the thermodynamics therefore advances the two copies on
     * alternate steps, each from where it was two steps earlier, and the
     * surface fluxes of one step are computed from the other copy. An initial
     * state read from a file goes into one copy only; the other starts
     * ice-free. True reproduces that, so a ROMS run without ICE_ADVECT can be
     * matched. False keeps one state, which is the physically sensible
     * thermodynamic model.
     */
    bool no_advect_two_states = true;

    //! Whether this configuration advects the ice
    bool advects () const { return do_advect; }

    /**
     * Ice and snow albedo: false is the default formulation of ROMS
     * ice_albedo.F, true is Ebert and Curry (1992), ROMS ICE_ALB_EC92. Set
     * by remora.ice.albedo_type = default | ec92. LAKE_ICE uses ec92.
     */
    bool albedo_ec92 = false;

    /**
     * Use the model's own rho0 and Cp wherever ROMS hard-codes a water
     * density or heat capacity in the ice module.
     *
     * ROMS fixes RhoSW = 1026 kg/m3, Cp_w = 3990 J/kg/K and
     * RhoCpr = 0.2442754e-6 = 1/(1026*3990) as PARAMETERS in ice_mk.h,
     * whatever rho0 the run is using, and writes `1000.0 + rho` in one
     * branch of ice_frazil.F against `rho0 + rho` in the other. On fresh
     * water with rho0 = 1000 that is:
     *
     *     RhoSW    1026     vs rho0 1000        +2.60 %
     *     Cp_w     3990     vs Cp   3985        +0.13 %
     *     RhoCpr   2.4428e-7 vs 1/(rho0*Cp)     -2.66 %
     *
     * RhoCpr multiplies every ice-derived surface heat flux.
     *
     * The default is false, which reproduces ROMS. Set it true for a run whose
     * ice uses the same water density and heat capacity as its ocean.
     */
    bool use_rho0 = false;

    /**
     * The ice variables that carry a lateral boundary condition, in the order
     * the bc_type table is indexed. ROMS applies ice_bc2d.F to the scalars,
     * ice_tibc.F to Tice and ice_uibc/ice_vibc.F to the velocities, so the
     * routine each one uses is implied by its entry here.
     */
    static const char** ice_bc_var_names (int& n);

    /**
     * Which row of bc_type/bc_value carries this IceS component's condition, or -1 for a
     * component that has none.
     */
    static int ice_bc_var_for_comp (int comp);

    /**
     * The ROMS varinfo.yaml name of the time coordinate for this variable's
     * boundary records, and the stagger of those records.
     *
     * varinfo.yaml gives `Aice_time` for the concentration and `ocean_time`
     * for every other ice boundary variable -- an inconsistency in ROMS, not
     * a design, but a file written for ROMS will follow it. Overridable per
     * variable with `remora.ice.bc_time_name.<var>`.
     */
    static const char* ice_bc_time_name (int bc_var);

    /** 0 = rho point, 1 = u point, 2 = v point, for the boundary records. */
    static int ice_bc_stagger (int bc_var);

    /**
     * Does this side of this variable need boundary DATA?
     *
     * ROMS's `%acquire`, set in inp_decode.F:1779-1805 for exactly three of
     * its LBC keywords: CLA (clamped), MIX (mixed) and RADNUD (radiation with
     * nudging). A plain RAD acquires nothing -- the radiation condition is
     * self-contained -- which is why nudging has to be asked for separately
     * before a file is required.
     *
     * In ROMS `nudging` is only ever set by RADNUD, so it cannot occur
     * without radiation. This port carries nudging as an independent per-side
     * switch, and the kernels only ever use it inside the radiation branch, so
     * the two agree: nudging without radiation changes nothing and needs no
     * data.
     */
    bool ice_bc_needs_data (int bc_var, int side) const
    {
        const IceBCType t = bc_type[bc_var][side];
        return (t == IceBCType::clamped) || (t == IceBCType::mixed)
            || ((t == IceBCType::radiation) && bc_nudging[side]);
    }

    //! [variable][side] lateral condition. Defaults to ROMS's `closed`,
    //! which is a gradient condition and is what a closed basin wants.
    IceBCType bc_type[16][4] = {{IceBCType::closed}};

    /**
     * The boundary VALUE each side holds, per variable, in ROMS's side order.
     *
     * `clamped`, `mixed` and the nudging term all need one; ROMS reads them
     * from BOUNDARY(ng)%ice_west/east/south/north, filled from a boundary
     * NetCDF file. REMORA reads such a file too (REMORA_ice_bdry_data.cpp).
     * This constant is used for a variable the file does not hold, or with no
     * file:
     *
     *   remora.ice.bc_value.Aice = 0.0 0.0 0.9 0.0
     *
     * One entry is broadcast to all four sides, and the default is zero. With every side
     * `closed` (the default) the values are never read.
     */
    amrex::Real bc_value[16][4] = {{amrex::Real(0.0)}};

    //! Per-side nudging, which in ROMS modifies `radiation` and is not a type.
    bool bc_nudging[4] = {false, false, false, false};
    //! Inflow and outflow nudging timescales [1/s] (ROMS Tobc_in/Tobc_out).
    amrex::Real bc_tau_in  = amrex::Real(0.0);
    amrex::Real bc_tau_out = amrex::Real(0.0);

    /**
     * ROMS's gamma2: the tangential-velocity slip coefficient at a `closed`
     * ice boundary. +1 free slip, -1 no slip. ROMS shares roms.in's value
     * between the ocean and the ice; this port keeps `bc_gamma2` separate so
     * the ice can be set without touching REMORA's own momentum boundary.
     */
    amrex::Real bc_gamma2 = amrex::Real(1.0);

    /**
     * The dt the nudging timescales are multiplied by (ROMS `tau*dt(ng)`).
     * Set by the driver to the ice step before the conditions are applied,
     * rather than threaded through every call site.
     */
    amrex::Real bc_dt = amrex::Real(0.0);

    //! Reference density [kg/m3], copied from remora.rho0 at init.
    amrex::Real rho0 = amrex::Real(1025.0);
    //! Sea-water heat capacity [J/kg/K], copied from REMORA's Cp at init.
    amrex::Real Cp = amrex::Real(3985.0);

    //! Water density the ice thermodynamics should use.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real RhoSW_eff () const {
        return use_rho0 ? rho0 : amrex::Real(1026.0);
    }
    //! Water heat capacity the ice thermodynamics should use.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real Cp_w_eff () const {
        return use_rho0 ? Cp : amrex::Real(3990.0);
    }
    //! 1/(rho*Cp), the W/m2 -> C m/s conversion the ice fluxes use.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real RhoCpr_eff () const {
        return use_rho0 ? (amrex::Real(1.0) / (rho0 * Cp))
                        : amrex::Real(0.2442754e-6);
    }

    //! Emissivity of the ice surface (ROMS ice_emiss).
    amrex::Real ice_emiss = amrex::Real(0.97);

    //! Specific heat of air [J/kg/K] (ROMS spec_heat_air).
    amrex::Real spec_heat_air = amrex::Real(1004.0);

    //! Air-ice heat transfer coefficient (ROMS trans_coeff).
    amrex::Real trans_coeff = amrex::Real(1.75e-3);

    //! Latent heat of sublimation [J/kg] (ROMS sublimation).
    amrex::Real sublimation = amrex::Real(2.834e6);

    /**
     * Use the evidently intended Celsius comparison in the ice albedo.
     *
     * false (default) reproduces ROMS, where the melt branch is unreachable
     * because a Celsius temperature is tested against a Kelvin threshold, so
     * the albedo is always dry -- 0.85 over snow, 0.65 over bare ice. See
     * ice_albedo() in REMORA_IceBulkFlux.H.
     */
    bool albedo_melt_fix = false;

    /**
     * Debug probe. When probe_i and probe_j are both >= 0, the driver prints
     * that cell's column every probe_int steps (a line starting ICEPROBE): the
     * inputs the model hands the thermodynamics, followed by the state it
     * produced. Several regression tests read their value from this line.
     *
     * Off by default and prints nothing when off.
     */
    int probe_i = -1;
    int probe_j = -1;
    int probe_int = 1;

    /**
     * Also print the GHOST values the lateral boundary conditions wrote, on
     * the four sides at the probe's own row and column, every probe_int
     * steps.
     *
     * The probe itself reports an interior cell, which the conditions do not write.
     */
    bool probe_bc = false;

    /**
     * @param[in ] remora_prefix the ParmParse prefix (`remora`), so the
     *             sea-ice keys read as `remora.ice.*`
     * @param[out] host_bc_time_name per-variable override of the boundary
     *             time-coordinate name (`remora.ice.bc_time_name.<var>`),
     *             written here but held by the caller, not by this struct --
     *             see the static_assert below.
     */
    void init_params (const std::string& remora_prefix,
                       std::string (&host_bc_time_name)[16]);
};

// SeaIceParameters is captured by value into device lambdas, so every member must be
// device-copyable. The time-coordinate names are strings, which are not, so they are a
// member of REMORA.
static_assert(std::is_trivially_copyable<SeaIceParameters>::value,
              "SeaIceParameters is captured into device kernels");

/**
 * The sea-ice fields offered to plotfile output, under the names ROMS's varinfo.yaml gives
 * them.
 *
 * Each entry is {name, which MultiFab, which component}: `in_state` true means
 * the prognostic Si, false the internal Fi. set2DPlotVariables() and the plotfile
 * writer both read this table, so the names they accept and write agree.
 */
struct IcePlotVar { const char* name; bool in_state; int comp; };

inline const IcePlotVar* ice_plot_vars (int& n)
{
    static const IcePlotVar v[] = {
        {"Aice",            true,  IceS::Aice},
        {"ice_thickness",   true,  IceS::Hice},
        {"snow_thickness",  true,  IceS::Hsno},
        {"meltpond_thickness", true, IceS::Hmel},
        {"ice_age",         true,  IceS::Iage},
        {"Tice",            true,  IceS::Tice},
        {"ice_sst",         false, IceF::Isst},
        {"under_ice_temp",  false, IceF::T0mk},
        {"under_ice_salt",  false, IceF::S0mk},
        {"io_massflux",     false, IceF::IOmf},
    };
    n = static_cast<int>(sizeof(v) / sizeof(v[0]));
    return v;
}

/**
 * The rho-point ice fields the NetCDF history file carries, under the names
 * ROMS varinfo.yaml gives them, so a ROMS history file and a REMORA one can
 * be compared variable by variable.
 *
 * The units are varinfo.yaml's with three exceptions, where ROMS's label does
 * not describe the number both models write. ice_age is labelled in seconds
 * there and in ice_mod.h, but ice_mk.h advances it by dtice*sec2day: days. The
 * thickness linked with it is then in meter days.
 * ice_divergence is labelled a divergence rate in m3/s, but ice_smolar.h
 * stores the thickness change the advection made, divided by the time step:
 * m/s, positive where the ice converges.
 */
struct IceHistVar { const char* name; const char* long_name; const char* units; bool in_state; int comp; };

inline const IceHistVar* ice_hist_vars (int& n)
{
    static const IceHistVar v[] = {
        {"Aice",               "fraction of cell covered by ice",       "nondimensional",   true,  IceS::Aice},
        {"ice_thickness",      "average ice thickness in cell",         "meter",            true,  IceS::Hice},
        {"meltpond_thickness", "melt pond water thickness on ice",      "meter",            true,  IceS::Hmel},
        {"snow_thickness",     "average thickness of snow coverage",    "meter",            true,  IceS::Hsno},
        {"ice_age",            "age of the ice",                        "day",              true,  IceS::Iage},
        {"ice_Sxx",            "internal ice stress, xx-component",     "newton meter-1",   true,  IceS::ISxx},
        {"ice_Sxy",            "internal ice stress, xy-component",     "newton meter-1",   true,  IceS::ISxy},
        {"ice_Syy",            "internal ice stress, yy-component",     "newton meter-1",   true,  IceS::ISyy},
        {"Tice",               "interior ice temperature",              "Celsius",          true,  IceS::Tice},
        {"ice_enthalpy",       "enthalpy of the ice/brine system",      "Celsius meter",    true,  IceS::Enth},
        {"ice_age_thickness",  "thickness associated with ice age",     "meter day",          true,  IceS::Hage},
        {"io_friction_vel",    "ice-ocean friction velocity",           "meter second-1",   false, IceF::IOfv},
        {"io_massflux",        "ice-ocean mass flux",                   "meter second-1",   false, IceF::IOmf},
        {"io_transfer_coef",   "ice-ocean momentum transfer coefficient","meter second-1",  false, IceF::IOmt},
        {"ice_sst",            "temperature at the ice/snow surface",   "Celsius",          false, IceF::Isst},
        {"ice_snow_heat_conductivity", "ice/snow heat conductivity",    "watt meter-2 kelvin-1", false, IceF::Qcon},
        {"ice_snow_heat_flux", "surface heat flux right-hand side over ice/snow", "watt meter-2", false, IceF::Qrhs},
        {"under_ice_salt",     "salinity of the molecular sublayer under ice", "nondimensional", false, IceF::S0mk},
        {"under_ice_temp",     "temperature of the molecular sublayer under ice", "Celsius", false, IceF::T0mk},
        {"ice_divergence",     "rate of change of ice thickness by advection", "meter second-1", false, IceF::Wdiv},
        {"ai_melt_freeze",     "rate of melt/freeze at the air-ice interface",   "meter second-1", false, IceF::W_ai},
        {"ao_melt_freeze",     "rate of melt/freeze at the air-ocean interface", "meter second-1", false, IceF::W_ao},
        {"frazil_freeze",      "rate of ice accretion by frazil growth",         "meter second-1", false, IceF::W_fr},
        {"io_melt_freeze",     "rate of melt/freeze at the ice-ocean interface", "meter second-1", false, IceF::W_io},
    };
    n = static_cast<int>(sizeof(v) / sizeof(v[0]));
    return v;
}

/**
 * Physical constants of the thermodynamics, from ROMS ice_mk.h. These are
 * parameters there, not inputs, so they are constants here too.
 *
 * Three of them collide with a constant REMORA already defines under a
 * different value, and the ice module deliberately keeps ROMS's:
 *
 *   - von Karman: ice_mk.h declares a local `kappa = 0.4`; REMORA's
 *     `vonKar` is 0.41. The ice sublayer formulation is Mellor-Kantha's,
 *     tuned with 0.4, so `IceKappa` below is 0.4.
 *   - water heat capacity: ice_mk.h uses `Cp_w = 3990`; REMORA's `Cp` is
 *     3985. `Cp_w` below is 3990.
 *   - specific heat of air: ice.in ships 1004.0; REMORA's `blk_Cpa` is
 *     1004.67 (Businger 1982). The ice value is a runtime parameter
 *     (`spec_heat_air`) defaulting to 1004.0.
 *
 */
namespace Thermo {
    constexpr amrex::Real AlphIc   = amrex::Real(2.034);        //!< ice conductivity [W/m/K]
    constexpr amrex::Real AlphSn   = amrex::Real(0.31);         //!< snow conductivity [W/m/K]
    constexpr amrex::Real Cp_i     = amrex::Real(2093.0);       //!< ice heat capacity [J/kg/K]
    // RhoSW, Cp_w and RhoCpr are DELIBERATELY not here: ROMS fixes them as
    // parameters regardless of rho0, and whether to honour that is a runtime
    // choice. Use SeaIceParameters::RhoSW_eff(), Cp_w_eff() and RhoCpr_eff().
    constexpr amrex::Real eps      = amrex::Real(1.0e-4);       //!< zero-division guard
    constexpr amrex::Real frln     = amrex::Real(-0.0543);      //!< freezing-line slope [C/psu]
    constexpr amrex::Real hfus     = amrex::Real(3.347e5);      //!< latent heat of fusion [J/kg]
    constexpr amrex::Real IceKappa = amrex::Real(0.4);          //!< von Karman, ROMS ice value; NOT REMORA's vonKar
    constexpr amrex::Real nu       = amrex::Real(1.8e-6);       //!< kinematic viscosity [m2/s]
    constexpr amrex::Real prs      = amrex::Real(2432.0);       //!< salt Schmidt number
    constexpr amrex::Real prt      = amrex::Real(13.0);         //!< temperature Prandtl number
    constexpr amrex::Real sice_ref = amrex::Real(3.2);          //!< reference ice salinity [psu]
    constexpr amrex::Real tpr      = amrex::Real(0.85);         //!< turbulent Prandtl number
    // The Yaglom-Kader coefficient. ROMS declares it `real(r8), parameter :: ykf = 3.14`, the
    // only parameter of its block without an `_r8` suffix, so its value is the single-precision
    // 3.14 promoted to double. This is that value.
    constexpr amrex::Real ykf      = amrex::Real(3.1400001049041748046875);
    constexpr amrex::Real z0ii     = amrex::Real(0.02);         //!< ice roughness length [m]
}

} // namespace REMORASeaIce

#endif
