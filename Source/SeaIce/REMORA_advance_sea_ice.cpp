#include <REMORA.H>
#include <REMORA_IceThermo.H>
#include <REMORA_IceAdvect.H>
#include <AMReX_FabArrayUtility.H>
#include <AMReX_ParallelDescriptor.H>

#include <sstream>
#include <string>
#include <vector>

using namespace amrex;

namespace {
/**
 * \brief A host-readable copy of a (possibly device) MultiFab's data.
 *
 * The ice probe prints MultiFab data from the host, which device memory does not allow, so
 * it reads a pinned copy.
 */
MultiFab
host_probe_copy (const MultiFab& src)
{
    MultiFab dst(src.boxArray(), src.DistributionMap(), src.nComp(),
                 src.nGrowVect(), MFInfo().SetArena(The_Pinned_Arena()));
    dtoh_memcpy(dst, src);
    Gpu::streamSynchronize();
    return dst;
}

/**
 * \brief Print the probe lines of all ranks from the IO rank, and empty the buffer.
 *
 * Ranks that print on their own can split each other's lines. Collective.
 */
void
print_probe_lines (std::ostringstream& lines)
{
    const std::string mine = lines.str();
    lines.str("");
    const int root = ParallelDescriptor::IOProcessorNumber();
    const std::vector<int> counts = ParallelDescriptor::Gather(static_cast<int>(mine.size()), root);
    std::vector<int> offsets(counts.size(), 0);
    for (std::size_t n = 1; n < counts.size(); ++n) { offsets[n] = offsets[n-1] + counts[n-1]; }
    std::vector<char> all(counts.empty() ? 1 : offsets.back() + counts.back() + 1);
    ParallelDescriptor::Gatherv(mine.data(), static_cast<int>(mine.size()), all.data(),
                                counts, offsets, root);
    if (ParallelDescriptor::IOProcessor()) {
        Print() << std::string(all.data(), all.size() - 1);
    }
}
}

/**
 * \brief Advance the sea-ice model one baroclinic step.
 *
 * Called from setup_step after the surface fluxes are assembled and before
 * the ocean right-hand sides, which is where ROMS main3d.F steps its ice:
 * bulk_flux, set_vbc, seaice, then rhs3d. The ice therefore reads the
 * start-of-step ocean (ROMS nrhs) and the surface stress and tracer fluxes it
 * writes drive the ocean in the same step.
 *
 * The stage order is ROMS seaice.F: ice_spdiw, ice_vbc, the EVP subcycle,
 * the advection, then the thermodynamics.
 */
void
REMORA::advance_sea_ice (int lev, Real time, Real dt_lev)
{
    BL_PROFILE("REMORA::advance_sea_ice()");
    if (!REMORASeaIce::has_sea_ice(sea_ice_model)) { return; }

    // Without the bulk flux there is no air-ice flux and no open-water flux to weight by the
    // ice cover, and the ice would replace the driver's fluxes with zero.
    if (running_with_coupling_driver && !DriverUsesStateForcing(driver_atmos_forcing_mode)) {
        Abort("remora.ice_model needs REMORA's bulk fluxes, but the coupling driver passes"
              " fluxes, so REMORA computes none");
    }

    // ROMS sets dtice to the baroclinic step
    const Real dtice = dt_lev;
    // What the boundary nudging timescales get multiplied by (ROMS's `tau*dt(ng)`)
    sea_ice_params.bc_dt = dtice;

    // Once, before any stage, as ROMS set_data: every condition in the step reads the same
    // boundary values.
    update_ice_bdry_data(lev, time);

    const auto& p = sea_ice_params;
    const int N = geom[lev].Domain().bigEnd(2);

    // The start-of-step ocean state, ROMS t(:,:,:,nrhs,:)
    MultiFab& S_old = *cons_old[lev];

    // ice_spdiw and ice_vbc are outside ROMS's ICE_MOMENTUM guard, so they run whether or not
    // the ice moves: ice_vbc replaces the wind stress under the ice.
    // With one ice state, keep the state the step starts from: it is the old
    // level the radiation and mixed boundary conditions read (ROMS know = liold)
    // in the advection and again after the thermodynamics.
    const bool one_state = p.advects() || !p.no_advect_two_states;
    if (one_state) {
        MultiFab::Copy(*vec_ice_s_old[lev], *vec_ice_s_new[lev], 0,0,
                       REMORASeaIce::IceS::NumComps, vec_ice_s_new[lev]->nGrowVect());
    }

    ice_spdiw(lev);
    ice_vbc(lev);
    // The EVP dynamics, which is a no-op unless remora.ice.do_momentum.
    ice_dynamics(lev, dtice);
    // ROMS seaice.F order: the enthalpy, the advection and its limiter come
    // after the dynamics and before the thermodynamics.
    ice_advect(lev, dtice);
    // The thermodynamics updates the current ice state in place (ROMS linew).
    // Without an advection ROMS has two states and has just switched which is current,
    // so the thermodynamics picks up the one it left two steps ago; the
    // internal stress is the exception, since ROMS indexes it separately.
    if (!one_state) {
        std::swap(vec_ice_s_old[lev], vec_ice_s_new[lev]);
        MultiFab::Copy(*vec_ice_s_new[lev], *vec_ice_s_old[lev], REMORASeaIce::IceS::ISxx,
                       REMORASeaIce::IceS::ISxx, 3, vec_ice_s_new[lev]->nGrowVect());
    }

    // The probe reads host copies, made before the MFIter loop because an MFIter cannot be
    // nested in another.
    const bool l_probe = (p.probe_i >= 0 && p.probe_j >= 0 &&
                          (istep[lev] % p.probe_int == 0));
    MultiFab h_s, h_f, h_cons, h_su, h_sv, h_qi, h_qo, h_sn, h_rn, h_zr, h_zw;
    std::ostringstream probe_lines;
    probe_lines.precision(17);
    if (l_probe) {
        h_s    = host_probe_copy(*vec_ice_s_new[lev]);
        h_f    = host_probe_copy(*vec_ice_f[lev]);
        h_cons = host_probe_copy(S_old);
        h_su   = host_probe_copy(*vec_sustr[lev]);
        h_sv   = host_probe_copy(*vec_svstr[lev]);
        h_qi   = host_probe_copy(*vec_Qnet_ai[lev]);
        h_qo   = host_probe_copy(*vec_Qnet_ao[lev]);
        h_sn   = host_probe_copy(*vec_snow[lev]);
        h_rn   = host_probe_copy(*vec_rain[lev]);
        h_zr   = host_probe_copy(*vec_z_r[lev]);
        h_zw   = host_probe_copy(*vec_z_w[lev]);
    }

    // ROMS ICE_THERMO. Without it the ocean keeps the surface fluxes the bulk
    // flux computed, and the ice state is whatever the dynamics left.
    if (p.do_thermo) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*vec_ice_s_new[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real      > const& Si     = vec_ice_s_new[lev]->array(mfi);
        Array4<Real      > const& Fi     = vec_ice_f[lev]->array(mfi);
        // The applied surface tracer flux, ROMS FORCES%stflx after set_vbc
        Array4<Real      > const& stflx  = vec_stflx[lev]->array(mfi);
        Array4<Real const> const& cons   = S_old.const_array(mfi);
        Array4<Real const> const& sustr  = vec_sustr[lev]->const_array(mfi);
        Array4<Real const> const& svstr  = vec_svstr[lev]->const_array(mfi);
        Array4<Real const> const& Qai    = vec_Qnet_ai[lev]->const_array(mfi);
        Array4<Real const> const& Qao    = vec_Qnet_ao[lev]->const_array(mfi);
        Array4<Real const> const& snow   = vec_snow[lev]->const_array(mfi);
        Array4<Real const> const& rain   = vec_rain[lev]->const_array(mfi);
        Array4<Real const> const& z_r    = vec_z_r[lev]->const_array(mfi);
        Array4<Real const> const& z_w    = vec_z_w[lev]->const_array(mfi);
        Array4<Real const> const& mskr   = vec_mskr[lev]->const_array(mfi);

        Box bx = mfi.tilebox();
        bx.setSmall(2,0); bx.setBig(2,0);

        // Debug probe: print one cell's thermodynamic inputs and result.
        // Host-side and outside the ParallelFor, so it costs nothing when off.
        // On the GPU the data is device memory, so it is read from the pinned
        // host copies made before this loop (host_probe_copy).
        const int pi = p.probe_i, pj = p.probe_j;
        if (l_probe && bx.contains(IntVect(pi,pj,0)))
        {
            const auto sa = h_s.const_array(mfi);
            const auto fa = h_f.const_array(mfi);
            const auto ca = h_cons.const_array(mfi);
            const auto su = h_su.const_array(mfi);
            const auto sv = h_sv.const_array(mfi);
            const auto qi = h_qi.const_array(mfi);
            const auto qo = h_qo.const_array(mfi);
            const auto sn = h_sn.const_array(mfi);
            const auto rn = h_rn.const_array(mfi);
            const auto zr = h_zr.const_array(mfi);
            const auto zw = h_zw.const_array(mfi);
            const Real utau = REMORASeaIce::ice_utau(su(pi,pj,0), su(pi+1,pj,0),
                                                     sv(pi,pj,0), sv(pi,pj+1,0));
            // Printed as two equal stress components whose utau is this one
            const Real s_eq = utau * utau / std::sqrt(Real(2.0));
            probe_lines
                << "ICEPROBE step " << istep[lev] << " cell (" << pi << "," << pj << ") IN "
                << ca(pi,pj,N,Temp_comp) << " " << ca(pi,pj,N,Salt_comp) << " "
                << s_eq << " " << s_eq << " "
                << qi(pi,pj,0) << " " << qo(pi,pj,0) << " "
                << sn(pi,pj,0) << " " << rn(pi,pj,0) << " "
                << (zw(pi,pj,N+1) - zr(pi,pj,N)) << " "
                << sa(pi,pj,0,REMORASeaIce::IceS::Aice) << " "
                << sa(pi,pj,0,REMORASeaIce::IceS::Hice) << " "
                << sa(pi,pj,0,REMORASeaIce::IceS::Hsno) << " "
                << sa(pi,pj,0,REMORASeaIce::IceS::Hmel) << " "
                << sa(pi,pj,0,REMORASeaIce::IceS::Iage) << " "
                << sa(pi,pj,0,REMORASeaIce::IceS::Tice) << " "
                << fa(pi,pj,0,REMORASeaIce::IceF::Qcon) << " "
                << fa(pi,pj,0,REMORASeaIce::IceF::Qrhs) << " "
                << fa(pi,pj,0,REMORASeaIce::IceF::T0mk) << " "
                << fa(pi,pj,0,REMORASeaIce::IceF::S0mk) << " "
                << fa(pi,pj,0,REMORASeaIce::IceF::W_fr) << " "
                << dtice << "\n";
        }

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            // Land cells hold uninitialised memory on the GPU, so the
            // mask gates the whole column rather than only the output.
            if (mskr(i,j,0) < Real(0.5)) { return; }

            REMORASeaIce::IceColumn c;
            c.temp_top = cons(i,j,N,Temp_comp);
            c.salt_top = cons(i,j,N,Salt_comp);
            // z_w(N+1) - z_r(N): the surface half-cell, which sets the
            // sublayer roughness ratio.
            c.dztop    = z_w(i,j,N+1) - z_r(i,j,N);
            c.Qnet_ai  = Qai(i,j,0);
            c.Qnet_ao  = Qao(i,j,0);
            c.snow     = snow(i,j,0);
            c.rain     = rain(i,j,0);
            c.utau_in  = REMORASeaIce::ice_utau(sustr(i,j,0), sustr(i+1,j,0),
                                                svstr(i,j,0), svstr(i,j+1,0));

            c.Aice = Si(i,j,0,REMORASeaIce::IceS::Aice);
            c.Hice = Si(i,j,0,REMORASeaIce::IceS::Hice);
            c.Hsno = Si(i,j,0,REMORASeaIce::IceS::Hsno);
            c.Hmel = Si(i,j,0,REMORASeaIce::IceS::Hmel);
            c.Iage = Si(i,j,0,REMORASeaIce::IceS::Iage);
            c.Tice = Si(i,j,0,REMORASeaIce::IceS::Tice);

            c.Isst = Fi(i,j,0,REMORASeaIce::IceF::Isst);
            c.T0mk = Fi(i,j,0,REMORASeaIce::IceF::T0mk);
            c.S0mk = Fi(i,j,0,REMORASeaIce::IceF::S0mk);
            c.Qcon = Fi(i,j,0,REMORASeaIce::IceF::Qcon);
            c.Qrhs = Fi(i,j,0,REMORASeaIce::IceF::Qrhs);
            c.W_fr = Fi(i,j,0,REMORASeaIce::IceF::W_fr);

            // ROMS ice_mk.h overwrites the temperature flux and adds to the
            // salt flux, which set_vbc has already multiplied by the surface
            // salinity.
            c.stflx_temp = stflx(i,j,0,Temp_comp);
            c.stflx_salt = stflx(i,j,0,Salt_comp);

            REMORASeaIce::ice_thermo_column(c, p, dtice, mskr(i,j,0));
            // seaice.F calls ice_limit straight after the thermodynamics
            REMORASeaIce::ice_limit_state(c.Aice, c.Hice, c.Hsno, c.Hmel, c.Tice, c.Iage, p);

            Si(i,j,0,REMORASeaIce::IceS::Aice) = c.Aice;
            Si(i,j,0,REMORASeaIce::IceS::Hice) = c.Hice;
            Si(i,j,0,REMORASeaIce::IceS::Hsno) = c.Hsno;
            Si(i,j,0,REMORASeaIce::IceS::Hmel) = c.Hmel;
            Si(i,j,0,REMORASeaIce::IceS::Iage) = c.Iage;
            Si(i,j,0,REMORASeaIce::IceS::Tice) = c.Tice;

            Fi(i,j,0,REMORASeaIce::IceF::Isst) = c.Isst;
            Fi(i,j,0,REMORASeaIce::IceF::T0mk) = c.T0mk;
            Fi(i,j,0,REMORASeaIce::IceF::S0mk) = c.S0mk;
            // ice_mk.h adds the conduction through the ice to both in place.
            // bulk_fluxes rebuilds them next step, so this only makes the
            // stored fields the ones ROMS stores.
            Fi(i,j,0,REMORASeaIce::IceF::Qcon) = c.Qcon;
            Fi(i,j,0,REMORASeaIce::IceF::Qrhs) = c.Qrhs;
            Fi(i,j,0,REMORASeaIce::IceF::IOmf) = c.IOmf;
            Fi(i,j,0,REMORASeaIce::IceF::W_ai) = c.W_ai;
            Fi(i,j,0,REMORASeaIce::IceF::W_ao) = c.W_ao;
            Fi(i,j,0,REMORASeaIce::IceF::W_io) = c.W_io;
            Fi(i,j,0,REMORASeaIce::IceF::W_ro) = c.W_ro;

            // The full ROMS surface heat flux, shortwave included. With
            // remora.shortwave_penetration the shortwave that passes the ice is
            // then carried down the column in prestep, as ROMS SOLAR_SOURCE does.
            stflx(i,j,0,Temp_comp) = c.stflx_temp;
            stflx(i,j,0,Salt_comp) = c.stflx_salt;
        });
    }
    if (l_probe) { print_probe_lines(probe_lines); }

    // ice_mk.h closes the thermodynamics with the lateral conditions
    for (int c : {REMORASeaIce::IceS::Aice, REMORASeaIce::IceS::Hice,
                  REMORASeaIce::IceS::Hsno, REMORASeaIce::IceS::Hmel,
                  REMORASeaIce::IceS::Iage}) {
        ice_bc2d(lev, c);
    }
    ice_tibc(lev);
    }

    // ROMS also fills the physical ghost cells of Isst, Qcon, Qrhs and the surface fluxes
    // (bc_r2d_tile). Nothing reads them, so they are not filled here.
    if (l_probe) {
        h_s = host_probe_copy(*vec_ice_s_new[lev]);
        h_f = host_probe_copy(*vec_ice_f[lev]);
        const MultiFab h_u = host_probe_copy(*vec_ice_u[lev]);
        const MultiFab h_v = host_probe_copy(*vec_ice_v[lev]);
        for (MFIter mfi(*vec_ice_s_new[lev]); mfi.isValid(); ++mfi) {
            if (!mfi.tilebox().contains(IntVect(p.probe_i,p.probe_j,0))) { continue; }
            const auto sa = h_s.const_array(mfi);
            const auto fa = h_f.const_array(mfi);
            const int pi = p.probe_i, pj = p.probe_j;
            probe_lines
                << "ICEPROBE step " << istep[lev] << " cell (" << pi << "," << pj << ") OUT"
                << " Aice " << sa(pi,pj,0,REMORASeaIce::IceS::Aice)
                << " Hice " << sa(pi,pj,0,REMORASeaIce::IceS::Hice)
                << " Hsno " << sa(pi,pj,0,REMORASeaIce::IceS::Hsno)
                << " Hmel " << sa(pi,pj,0,REMORASeaIce::IceS::Hmel)
                << " Iage " << sa(pi,pj,0,REMORASeaIce::IceS::Iage)
                << " Tice " << sa(pi,pj,0,REMORASeaIce::IceS::Tice)
                << " Isst " << fa(pi,pj,0,REMORASeaIce::IceF::Isst)
                << " T0mk " << fa(pi,pj,0,REMORASeaIce::IceF::T0mk)
                << " S0mk " << fa(pi,pj,0,REMORASeaIce::IceF::S0mk)
                << "\n";
            // The ice velocity, so the dynamics can be seen to be doing
            // something (it is identically zero without remora.ice.do_momentum).
            const auto uia = h_u.const_array(mfi);
            const auto via = h_v.const_array(mfi);
            probe_lines
                << "ICEVEL   step " << istep[lev] << " cell (" << pi << "," << pj << ")"
                << " Uice " << uia(pi,pj,0)
                << " Vice " << via(pi,pj,0)
                << " IOvs " << fa(pi,pj,0,REMORASeaIce::IceF::IOvs)
                << " IOmt " << fa(pi,pj,0,REMORASeaIce::IceF::IOmt)
                << "\n";
        }
    }

    // The ghost values the lateral conditions wrote. A periodic side holds whatever the last
    // exchange left there.
    if (l_probe && p.probe_bc)
    {
        const Box& dom = geom[lev].Domain();
        const int pi = p.probe_i, pj = p.probe_j;
        const int wi = dom.smallEnd(0) - 1, ei = dom.bigEnd(0) + 1;
        const int sj = dom.smallEnd(1) - 1, nj = dom.bigEnd(1) + 1;
        const IntVect cells[4] = {IntVect(wi,pj,0), IntVect(pi,sj,0),
                                  IntVect(ei,pj,0), IntVect(pi,nj,0)};
        const char* tag[4] = {"W", "S", "E", "N"};
        const MultiFab h_u = host_probe_copy(*vec_ice_u[lev]);
        const MultiFab h_v = host_probe_copy(*vec_ice_v[lev]);
        for (MFIter mfi(*vec_ice_s_new[lev]); mfi.isValid(); ++mfi) {
            const Box& vbx = mfi.validbox();
            const auto sa  = h_s.const_array(mfi);
            const auto uia = h_u.const_array(mfi);
            const auto via = h_v.const_array(mfi);
            for (int k = 0; k < 4; ++k) {
                // The box that holds the interior cell beside the ghost cell
                const IntVect inside = amrex::max(amrex::min(cells[k], dom.bigEnd()), dom.smallEnd());
                if (!vbx.contains(IntVect(inside[0],inside[1],0))) { continue; }
                const int gi = cells[k][0], gj = cells[k][1];
                // The velocity BC sets a FACE on the side where the component
                // is normal, so u is reported at the boundary face on W/E and
                // at the ghost row on S/N; v the other way round.
                const int ui_i = (k == 0) ? gi + 1 : gi;
                const int vi_j = (k == 1) ? gj + 1 : gj;
                probe_lines
                    << "ICEBC    step " << istep[lev] << " side " << tag[k]
                    << " cell (" << gi << "," << gj << ")"
                    << " Aice " << sa(gi,gj,0,REMORASeaIce::IceS::Aice)
                    << " Hice " << sa(gi,gj,0,REMORASeaIce::IceS::Hice)
                    << " Tice " << sa(gi,gj,0,REMORASeaIce::IceS::Tice)
                    << " Enth " << sa(gi,gj,0,REMORASeaIce::IceS::Enth)
                    << " Uice " << uia(ui_i,gj,0)
                    << " Vice " << via(gi,vi_j,0)
                    << "\n";
            }
        }
    }
    if (l_probe) { print_probe_lines(probe_lines); }

    vec_ice_s_new[lev]->FillBoundary(geom[lev].periodicity());
    vec_ice_f[lev]->FillBoundary(geom[lev].periodicity());
    vec_stflx[lev]->FillBoundary(geom[lev].periodicity());
}
