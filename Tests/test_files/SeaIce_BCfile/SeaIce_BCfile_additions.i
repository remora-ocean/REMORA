
# SeaIce_BCfile: SeaIce_BC with the ice boundary values read from a NetCDF file.
#
# seaice_bry_we.nc holds three records, at 0, 0.02 and 1 days. In each, the value at index n
# along the side (n = 0 is the ghost row for rho and u points, the first face for v points) is
#
#                     day 0    day 0.02   day 1      per index
#   Aice              0.20     0.30       0.80       0.011
#   ice_thickness     0.10     0.25       0.50       0.013
#   Tice             -1.0     -2.5       -6.0       -0.017
#   Uice              0.05     0.09       0.21       0.019
#   Vice             -0.04    -0.07      -0.15       0.023
#
# plus 0.02 on the west. Step 75 is day 0.0434, past the kink at day 0.02.
#
# The west Aice side is radiation with nudging, which needs a value; the other west sides need
# none, and the file has none for them.
remora.ice.nc_bdry_file = seaice_bry_we.nc

#                                  WEST      SOUTH    EAST    NORTH
remora.ice.bc.Aice               = radiation periodic clamped periodic
remora.ice.bc.ice_thickness      = closed    periodic clamped periodic
remora.ice.bc.snow_thickness     = closed    periodic closed  periodic
remora.ice.bc.meltpond_thickness = closed    periodic closed  periodic
remora.ice.bc.ice_age            = closed    periodic closed  periodic
remora.ice.bc.Tice               = closed    periodic clamped periodic
remora.ice.bc.Uice               = closed    periodic clamped periodic
remora.ice.bc.Vice               = closed    periodic clamped periodic
remora.ice.bc_nudging            = 1 0 0 0

# Zero, so that a variable that does not read the file is seen.
remora.ice.bc_value.Aice          = 0.0 0.0 0.0 0.0
remora.ice.bc_value.ice_thickness = 0.0 0.0 0.0 0.0
remora.ice.bc_value.Tice          = 0.0 0.0 0.0 0.0
