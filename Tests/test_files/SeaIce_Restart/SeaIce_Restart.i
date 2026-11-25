# Sea-ice restart.
#
# The SeaIce_BC case, 40 steps with a checkpoint at 20, run by add_test_restart_exact:
# straight through in runA, restarted from runA/chk00020 in runB, and the two
# plt00040 must agree. 40 steps are enough for ice to grow from open water, the
# EVP subcycle to run and the advection to carry ice across the boundary, the
# three stages that write ice state another stage reads on the next step.
#
# The comparison sees the ice through the ocean fields: without the ice in the checkpoint,
# temp differs by 2.3e-05.
#
# The surface is cold because the BoundaryLayer profile puts its warm end at the
# top (T = T0 + 6.25*val1), so T0 is pushed down to get a freezing surface. The
# deep water is unphysically cold.

# ------------------  INPUTS TO MAIN PROGRAM  -------------------
remora.prob_name = BoundaryLayer

max_step = 40

remora.omp_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_lo     =      0.     0.    -47.5
geometry.prob_hi     =  100000. 5000.       0.

amr.n_cell           =  39     4      30

remora.is_periodic = 0 1 0

remora.boundary_per_variable = 1
##                      WEST              SOUTH             EAST              NORTH
remora.bc.temp.type   =  slipwall         periodic          orlanski_rad      periodic
remora.bc.salt.type   =  slipwall         periodic          orlanski_rad      periodic
remora.bc.u.type      =  slipwall         periodic          orlanski_rad      periodic
remora.bc.v.type      =  slipwall         periodic          orlanski_rad      periodic
remora.bc.ubar.type   =  slipwall         periodic          outflow           periodic
remora.bc.vbar.type   =  slipwall         periodic          outflow           periodic
remora.bc.zeta.type   =  slipwall         periodic          outflow           periodic
remora.bc.tke.type    =  slipwall         periodic          outflow           periodic

# TIME STEP CONTROL
remora.fixed_dt       = 50.0 # Timestep size (seconds)
remora.ndtfast = 30

# DIAGNOSTICS & VERBOSITY
remora.sum_interval  = 1       # timesteps between integrated/max quantities, if remora.v > 0
remora.v             = 0       # verbosity in REMORA.cpp (0: none, 1: integrated quantities, etc, 2: print boxes)

remora.bottom_stress_type = "quadratic"

# CHECKPOINT FILES
remora.check_file      = chk        # root name of checkpoint file
remora.check_int       = 20          # the restart point this case exists for

# PLOTFILES
remora.plot_file     = plt        # prefix of plotfile name
remora.plot_int      = 40        # only the final one, which is what is compared
remora.plot_vars_3d  = salt temp x_velocity y_velocity
remora.plotfile_type = amrex

# SOLVER CHOICE
remora.tracer_horizontal_advection_scheme = "upstream3" # upstream3 or centered4

remora.horizontal_mixing_type = "constant"
remora.visc2 = 0.0

remora.grid_scale_type = "analytic"

remora.eos_type = "nonlinear"

remora.bulk_fluxes = true

remora.air_temperature = -15.0
remora.air_humidity = 0.776
remora.air_pressure = 1013.48
remora.cloud = 0.5
remora.rain = 0.002
remora.eminusp = true
remora.eminusp_correct_ssh = true

# Linear EOS parameters
remora.R0    = 1027.0  # background density value (Kg/m3) used in Linear Equation of State
remora.S0    = 1.0     # nearly fresh
remora.T0    = -6.60   # puts the SURFACE near -0.37 C, just below freezing
remora.rho0  = 1000.0  # freshwater

remora.theta_s = 5.0
remora.theta_b = 1.5
remora.tcline = 20.0

remora.vertical_mixing_type="gls"
remora.Akt_bak = 1.0e-6
remora.Akv_bak = 1.0e-5

# Coriolis params
remora.use_coriolis = true
remora.coriolis_type = beta_plane
remora.coriolis_f0 = 9.25e-4
remora.coriolis_beta = 0.0

# As in BoundaryLayer.i; the expected values were made with the cap off.
remora.gls_length_cap = false

# ---- sea ice ----
remora.ice_model          = mellor_kantha

# Frazil off: the whole column is below freezing, and the frazil would turn its supercooling
# into metres of ice in a step.
remora.ice.do_frazil      = false
remora.shortwave_penetration = true
remora.jerlov_type           = IB
remora.surface_radiation_flux = 18.8     # net shortwave, W/m2, a winter day

remora.plot_vars_2d = zeta Aice ice_thickness snow_thickness meltpond_thickness \
                      ice_age Tice ice_sst under_ice_temp under_ice_salt \
                      io_massflux

# ---- the ice lateral boundary ----
#
# One entry per side in ROMS's order: west south east north.
##                                     WEST       SOUTH      EAST       NORTH
remora.ice.bc.Aice          = radiation periodic clamped periodic
remora.ice.bc.ice_thickness = radiation periodic clamped periodic
remora.ice.bc.snow_thickness= radiation periodic clamped periodic
remora.ice.bc.meltpond_thickness = radiation periodic clamped periodic
remora.ice.bc.ice_age       = radiation periodic clamped periodic
remora.ice.bc.Tice          = gradient  periodic clamped periodic
remora.ice.bc.Uice          = radiation periodic closed  periodic
remora.ice.bc.Vice          = closed    periodic closed  periodic

# Ice coming IN from the east: half-covered, 20 cm thick, -4 C. The west is
# radiation, so its value is only the nudging target.
remora.ice.bc_value.Aice          = 0.0 0.0 0.5  0.0
remora.ice.bc_value.ice_thickness = 0.0 0.0 0.2  0.0
remora.ice.bc_value.Tice          = 0.0 0.0 -4.0 0.0

# Nudging on the radiation side, with ROMS's inflow/outflow pair of rates.
remora.ice.bc_nudging  = 1 0 0 0
remora.ice.bc_tau_in   = 1.0e-4
remora.ice.bc_tau_out  = 1.0e-5

# Free slip for the tangential ice velocity at a closed side (ROMS gamma2).
remora.ice.bc_gamma2 = 1.0
