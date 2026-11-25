#include <REMORA.H>
#include <REMORA_SeaIce.H>
#include <REMORA_IceBC.H>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

using namespace amrex;

namespace REMORASeaIce {

SeaIceModel
parse_sea_ice_model (const std::string& name)
{
    if (name == "none") { return SeaIceModel::none; }
    if (name == "mellor_kantha") { return SeaIceModel::mellor_kantha; }
    amrex::Abort("remora.ice_model must be none or mellor_kantha; got '" + name + "'");
    return SeaIceModel::none;
}

std::string
sea_ice_model_name (SeaIceModel m)
{
    switch (m) {
        case SeaIceModel::mellor_kantha: return "mellor_kantha";
        default:                         return "none";
    }
}

/**
 * The ice variables that carry a lateral boundary condition.
 *
 * The order fixes how bc_type is indexed. The names are the ROMS output names
 * already used for the plotfile, so a deck says the same thing in both
 * places. The routine each one uses in ROMS is implied: the scalars go
 * through ice_bc2d.F, Tice through ice_tibc.F, and the velocities through
 * ice_uibc.F / ice_vibc.F.
 */
const char**
SeaIceParameters::ice_bc_var_names (int& n)
{
    // ROMS's own names for its 15 ice state variables, in ROMS's own order,
    // from the table in ROMS/External/ice.in:270-287. Transcribing an ice.in
    // should be mechanical, so the names and the order are ROMS's.
    //
    // isUevp and isVevp are deliberately absent: they are the EVP subcycle's
    // own working velocities, local to ice_dynamics, and no boundary condition
    // is applied to them. ROMS lists them in ice.in and then never reads an
    // LBC record for them either.
    static const char* v[] = {
        "Aice",               // isAice  ice concentration
        "ice_thickness",      // isHice
        "meltpond_thickness", // isHmel  meltwater thickness on ice
        "snow_thickness",     // isHsno
        "ice_age",            // isIage
        "ice_Sxx",            // isISxx  internal ice xx-stress
        "ice_Sxy",            // isISxy
        "ice_Syy",            // isISyy
        "Tice",               // isTice  ice interior temperature
        "Uice",               // isUice
        "Vice",               // isVice
        "ice_enthalpy",       // isEnth
        "ice_age_thickness"   // isHage
    };
    n = static_cast<int>(sizeof(v) / sizeof(v[0]));
    return v;
}

const char*
SeaIceParameters::ice_bc_time_name (int bc_var)
{
    // varinfo.yaml, the `time:` field of the *_west/_east/_south/_north
    // entries. Aice is the odd one out and that is ROMS's, not ours.
    return (bc_var == 0) ? "Aice_time" : "ocean_time";
}

int
SeaIceParameters::ice_bc_stagger (int bc_var)
{
    // varinfo.yaml types every ice boundary record r2dvar except Uice and
    // Vice. It calls BOTH of those u2dvar (ROMS/External/varinfo.yaml:1953ff),
    // which for Vice is wrong -- a v-point record on an east or west boundary
    // has one more element than a u-point one. ROMS does not read the record
    // through that type (set_data.F:1251 passes explicit bounds, and uses
    // iceIstrM = 1 for isVice against 0 for the rest, which is the v stagger),
    // so the mistyping does not bite there; it would bite anything that did
    // trust the type. This port follows the bounds, not the type.
    int n = 0;
    const auto* v = ice_bc_var_names(n);
    if (bc_var < 0 || bc_var >= n) { return 0; }
    const std::string nm(v[bc_var]);
    if (nm == "Uice") { return 1; }
    if (nm == "Vice") { return 2; }
    return 0;
}

int
SeaIceParameters::ice_bc_var_for_comp (int comp)
{
    // Must stay in step with ice_bc_var_names above.
    switch (comp) {
        case IceS::Aice: return 0;    // "Aice"
        case IceS::Hice: return 1;    // "ice_thickness"
        case IceS::Hmel: return 2;    // "meltpond_thickness"
        case IceS::Hsno: return 3;    // "snow_thickness"
        case IceS::Iage: return 4;    // "ice_age"
        case IceS::ISxx: return 5;    // "ice_Sxx"
        case IceS::ISxy: return 6;    // "ice_Sxy"
        case IceS::ISyy: return 7;    // "ice_Syy"
        case IceS::Tice: return 8;    // "Tice"
        case IceS::Uice: return 9;    // "Uice"
        case IceS::Vice: return 10;   // "Vice"
        case IceS::Enth: return 11;   // "ice_enthalpy"
        case IceS::Hage: return 12;   // "ice_age_thickness"
        default: return -1;
    }
}

IceBCType
parse_ice_bc_type (const std::string& name)
{
    if (name == "radiation") return IceBCType::radiation;
    if (name == "clamped")   return IceBCType::clamped;
    if (name == "mixed")     return IceBCType::mixed;
    if (name == "closed")    return IceBCType::closed;
    if (name == "gradient")  return IceBCType::gradient;
    if (name == "periodic")  return IceBCType::periodic;
    amrex::Abort("ice boundary type must be one of closed, gradient, "
                 "radiation, clamped, mixed, periodic; got '" + name + "'");
    return IceBCType::closed;
}

std::string
ice_bc_type_name (IceBCType t)
{
    switch (t) {
        case IceBCType::radiation: return "radiation";
        case IceBCType::clamped:   return "clamped";
        case IceBCType::mixed:     return "mixed";
        case IceBCType::gradient:  return "gradient";
        case IceBCType::periodic:  return "periodic";
        default:                   return "closed";
    }
}

void
SeaIceParameters::init_params (const std::string& remora_prefix,
                                std::string (&host_bc_time_name)[16])
{
    // The names after the prefix are those of ROMS's ice.in.
    const std::string prefix = remora_prefix + ".ice";
    ParmParse pp(prefix);

    // ROMS reads these six from ice.in and, in the configuration ported here,
    // uses none of them: the air-ice fluxes take their constants from the bulk
    // flux, and no turning angle is applied. They are read so that an ice.in
    // can be transcribed, and a deck that sets one is told it has no effect.
    // Before the queryAdds below, which put every name in the table.
    {
        std::string unused;
        for (const char* name : {"AirRho", "stressAng", "ice_emiss", "spec_heat_air",
                                 "trans_coeff", "sublimation"}) {
            if (pp.contains(name)) { unused += std::string(" ") + name; }
        }
        if (!unused.empty() && ParallelDescriptor::IOProcessor()) {
            amrex::Warning("ice parameters with no effect, in REMORA as in ROMS:" + unused);
        }
    }

    pp.queryAdd("nEVP", nEVP);

    pp.queryAdd("AirRho", AirRho);
    pp.queryAdd("IceRho", IceRho);
    pp.queryAdd("SnowDryRho", SnowDryRho);
    pp.queryAdd("SnowWetRho", SnowWetRho);

    pp.queryAdd("Cd_ai", Cd_ai);
    pp.queryAdd("Cd_io", Cd_io);

    pp.queryAdd("Astrength", Astrength);
    pp.queryAdd("Pstar", Pstar);

    pp.queryAdd("min_ai", min_ai);
    pp.queryAdd("min_hi", min_hi);
    pp.queryAdd("max_ai", max_ai);
    pp.queryAdd("max_hmelt", max_hmelt);

    pp.queryAdd("zetaMin", zetaMin);
    pp.queryAdd("zetaMax", zetaMax);
    pp.queryAdd("stressAng", stressAng);
    pp.queryAdd("ellip_sq", ellip_sq);

    pp.queryAdd("ice_emiss", ice_emiss);
    pp.queryAdd("spec_heat_air", spec_heat_air);
    pp.queryAdd("trans_coeff", trans_coeff);
    pp.queryAdd("sublimation", sublimation);

    pp.queryAdd("albedo_melt_fix", albedo_melt_fix);
    pp.queryAdd("zmax_clamp", zmax_clamp);
    pp.queryAdd("pressure_freezing", pressure_freezing);
    // The ROMS-self-inconsistency switches. `bc_fix_all` sets every one, and
    // an individual key still wins over it, so `bc_fix_all = 1` with one
    // member set back to 0 is a usable combination.
    {
        bool all = false;
        pp.queryAdd("bc_fix_all", all);
        if (all) {
            bc_fix.tibc_mixed = true;
            bc_fix.uibc_north_nudging = true;
            bc_fix.vibc_west_closed = true;
            bc_fix.vibc_north_mask = true;
        }
        amrex::ParmParse ppf(prefix + ".bc_fix");
        ppf.queryAdd("tibc_mixed", bc_fix.tibc_mixed);
        ppf.queryAdd("uibc_north_nudging", bc_fix.uibc_north_nudging);
        ppf.queryAdd("vibc_west_closed", bc_fix.vibc_west_closed);
        ppf.queryAdd("vibc_north_mask", bc_fix.vibc_north_mask);
    }
    pp.queryAdd("do_frazil", do_frazil);
    pp.queryAdd("qnet_ao_open_water_sw", qnet_ao_open_water_sw);
    pp.queryAdd("mom_bulk", mom_bulk);
    pp.queryAdd("convsnow", convsnow);
    pp.queryAdd("do_thermo", do_thermo);
    pp.queryAdd("do_momentum", do_momentum);
    pp.queryAdd("advect_upwind", advect_upwind);
    pp.queryAdd("strength_quad", strength_quad);
    pp.queryAdd("do_advect", do_advect);
    pp.queryAdd("no_advect_two_states", no_advect_two_states);
    {
        std::string albedo_type = "default";
        pp.queryAdd("albedo_type", albedo_type);
        if (albedo_type == "ec92") {
            albedo_ec92 = true;
        } else if (albedo_type != "default") {
            // ROMS also has ALBEDO_CSIM, which reads a variable it never
            // computes and does not compile.
            amrex::Abort(prefix + ".albedo_type must be default or ec92; got '" + albedo_type + "'");
        }
    }
    pp.queryAdd("use_rho0", use_rho0);

    {
        int nv = 0;
        ice_bc_var_names(nv);
        for (int v = 0; v < nv; ++v) {
            for (int sd = 0; sd < IceSide::count; ++sd) {
                bc_type[v][sd]  = IceBCType::closed;
                bc_value[v][sd] = amrex::Real(0.0);
            }
        }
    }

    // Lateral boundary conditions, one entry per side in ROMS's own order --
    // west, south, east, north -- matching the LBC records in ice.in. A single
    // entry is broadcast to all four sides, which is how a uniformly closed
    // basin is written.
    //
    //   remora.ice.bc.Aice = closed closed closed closed
    //   remora.ice.bc.Uice = radiation closed radiation closed
    //
    // The default is `closed` everywhere, which is ROMS's gradient condition.
    {
        int n = 0;
        const auto* vars = ice_bc_var_names(n);
        // One condition for every ice variable; a bc.<var> below overrides it.
        if (pp.contains("bc_type")) {
            const int cnt = pp.countval("bc_type");
            if (cnt != 1 && cnt != IceSide::count) {
                amrex::Abort(prefix + ".bc_type needs 1 or 4 entries (west south east north)");
            }
            for (int sd = 0; sd < IceSide::count; ++sd) {
                std::string nm;
                pp.get("bc_type", nm, (cnt == 1) ? 0 : sd);
                const IceBCType t = parse_ice_bc_type(nm);
                for (int v = 0; v < n; ++v) { bc_type[v][sd] = t; }
            }
        }
        for (int v = 0; v < n; ++v) {
            const std::string key = std::string("bc.") + vars[v];
            if (!pp.contains(key.c_str())) { continue; }
            const int cnt = pp.countval(key.c_str());
            if (cnt != 1 && cnt != IceSide::count) {
                amrex::Abort(prefix + "." + key +
                             " needs 1 or 4 entries (west south east north)");
            }
            for (int sd = 0; sd < IceSide::count; ++sd) {
                std::string nm;
                pp.get(key.c_str(), nm, (cnt == 1) ? 0 : sd);
                bc_type[v][sd] = parse_ice_bc_type(nm);
            }
        }
        // The combinations ROMS's routines have no branch for. A velocity without a condition
        // leaves the boundary face to the momentum step, so it is refused; a Tice without one
        // keeps its ghost value.
        for (int v = 0; v < n; ++v) {
            const std::string var(vars[v]);
            const std::string key = prefix + ".bc." + var;
            for (int sd = 0; sd < IceSide::count; ++sd) {
                const IceBCType t = bc_type[v][sd];
                if ((var == "Uice" || var == "Vice") && t == IceBCType::mixed) {
                    amrex::Abort(key + ": mixed applies no condition to " + var
                                 + ", in REMORA as in ROMS");
                }
                if (var == "Tice" && t == IceBCType::radiation) {
                    if (ParallelDescriptor::IOProcessor()) {
                        amrex::Warning(key + ": radiation applies no condition to Tice,"
                                       " in REMORA as in ROMS");
                    }
                    break;
                }
            }
            if ((var == "ice_enthalpy" || var == "ice_age_thickness") &&
                pp.contains((std::string("bc.") + var).c_str()) &&
                ParallelDescriptor::IOProcessor()) {
                amrex::Warning(key + " has no effect: it takes its boundary values from "
                               + (var == "ice_enthalpy" ? "Tice" : "ice_age"));
            }
        }
        const std::string nkey = "bc_nudging";
        if (pp.contains(nkey.c_str())) {
            const int cnt = pp.countval(nkey.c_str());
            if (cnt != 1 && cnt != IceSide::count) {
                amrex::Abort(prefix + ".bc_nudging needs 1 or 4 entries");
            }
            for (int sd = 0; sd < IceSide::count; ++sd) {
                int b = 0;
                pp.get(nkey.c_str(), b, (cnt == 1) ? 0 : sd);
                bc_nudging[sd] = (b != 0);
            }
        }
        // The constant boundary values, same shape as bc.<var>.
        for (int v = 0; v < n; ++v) {
            const std::string key = std::string("bc_value.") + vars[v];
            if (!pp.contains(key.c_str())) { continue; }
            const int cnt = pp.countval(key.c_str());
            if (cnt != 1 && cnt != IceSide::count) {
                amrex::Abort(prefix + "." + key +
                             " needs 1 or 4 entries (west south east north)");
            }
            for (int sd = 0; sd < IceSide::count; ++sd) {
                amrex::Real x = amrex::Real(0.0);
                pp.get(key.c_str(), x, (cnt == 1) ? 0 : sd);
                bc_value[v][sd] = x;
            }
        }
        pp.queryAdd("bc_tau_in",  bc_tau_in);
        pp.queryAdd("bc_tau_out", bc_tau_out);
        pp.queryAdd("bc_gamma2", bc_gamma2);
        // Per-variable time-coordinate override. Host-only storage: see the
        // static_assert beside the struct definition.
        for (int v = 0; v < n; ++v) {
            const std::string key = std::string("bc_time_name.") + vars[v];
            pp.queryAdd(key.c_str(), host_bc_time_name[v]);
        }
    }

    pp.queryAdd("probe_i", probe_i);
    pp.queryAdd("probe_j", probe_j);
    pp.queryAdd("probe_int", probe_int);
    pp.queryAdd("probe_bc", probe_bc);
    if (probe_int < 1) { amrex::Abort(prefix + ".probe_int must be >= 1"); }

    if (min_ai > max_ai) {
        amrex::Abort(prefix + ".min_ai exceeds max_ai");
    }
    if (nEVP < 1) {
        amrex::Abort(prefix + ".nEVP must be at least 1");
    }
}

} // namespace REMORASeaIce

/**
 * Every ice MultiFab that has to survive a restart.
 *
 * The internal fields are included: some are written by one stage and read by an earlier
 * stage of the next step, among them the surface temperature, the molecular-sublayer state
 * and the ice-ocean transfer coefficient.
 *
 * Both levels of the state are here because the radiation and mixed conditions read the old
 * one (ROMS's `know`).
 */
amrex::Vector<std::pair<std::string, amrex::MultiFab*>>
REMORA::ice_checkpoint_fields (int lev)
{
    amrex::Vector<std::pair<std::string, amrex::MultiFab*>> v;
    if (!REMORASeaIce::has_sea_ice(sea_ice_model)) { return v; }

    // The prognostic ice state, both time levels.
    v.emplace_back("ice_s_old",  vec_ice_s_old[lev].get());
    v.emplace_back("ice_s_new",  vec_ice_s_new[lev].get());
    // The diagnostic/internal set: ROMS's Fi.
    v.emplace_back("ice_f",      vec_ice_f[lev].get());
    // The ice velocity, and the EVP subcycle's own pair.
    v.emplace_back("ice_u",      vec_ice_u[lev].get());
    v.emplace_back("ice_v",      vec_ice_v[lev].get());
    v.emplace_back("ice_u_prev", vec_ice_u_prev[lev].get());
    v.emplace_back("ice_v_prev", vec_ice_v_prev[lev].get());
    v.emplace_back("ice_uevp",   vec_ice_uevp[lev].get());
    v.emplace_back("ice_vevp",   vec_ice_vevp[lev].get());
    v.emplace_back("ice_uavg",   vec_ice_uavg[lev].get());
    v.emplace_back("ice_vavg",   vec_ice_vavg[lev].get());
    v.emplace_back("ice_aius",   vec_ice_aius[lev].get());
    v.emplace_back("ice_aivs",   vec_ice_aivs[lev].get());
    // Air-ice and air-ocean heat, the snowfall rate, the ice-attenuated shortwave, and the
    // open-water and air-ice stresses of the bulk flux.
    v.emplace_back("Qnet_ao",    vec_Qnet_ao[lev].get());
    v.emplace_back("Qnet_ai",    vec_Qnet_ai[lev].get());
    v.emplace_back("ice_snow",   vec_snow[lev].get());
    v.emplace_back("srflx_ice",  vec_srflx_ice[lev].get());
    v.emplace_back("sustr_ao",   vec_sustr_ao[lev].get());
    v.emplace_back("svstr_ao",   vec_svstr_ao[lev].get());
    v.emplace_back("sustr_ai",   vec_sustr_ai[lev].get());
    v.emplace_back("svstr_ai",   vec_svstr_ai[lev].get());
    // The ocean density anomaly, which the frazil reads
    if (sea_ice_params.do_frazil) { v.emplace_back("rho_ice", vec_rho_ice[lev].get()); }
    return v;
}
