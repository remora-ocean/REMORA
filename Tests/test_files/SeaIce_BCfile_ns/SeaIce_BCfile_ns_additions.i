
# SeaIce_BCfile_ns: SeaIce_BCfile with the south and north sides open instead of the west and
# east. The records run along x here, and the u record is the short one.
#
# seaice_bry_sn.nc holds the values listed in SeaIce_BCfile_additions.i, plus 0.03 on the south.
#
# The south Vice side is radiation with nudging, which needs a value; Aice is closed there and
# needs none.
remora.ice.nc_bdry_file = seaice_bry_sn.nc

remora.is_periodic = 0 0 0

#                       WEST      SOUTH     EAST          NORTH
remora.bc.temp.type   = slipwall  slipwall  orlanski_rad  slipwall
remora.bc.salt.type   = slipwall  slipwall  orlanski_rad  slipwall
remora.bc.u.type      = slipwall  slipwall  orlanski_rad  slipwall
remora.bc.v.type      = slipwall  slipwall  orlanski_rad  slipwall
remora.bc.ubar.type   = slipwall  slipwall  outflow       slipwall
remora.bc.vbar.type   = slipwall  slipwall  outflow       slipwall
remora.bc.zeta.type   = slipwall  slipwall  outflow       slipwall
remora.bc.tke.type    = slipwall  slipwall  outflow       slipwall

#                                  WEST    SOUTH     EAST    NORTH
remora.ice.bc.Aice               = closed  closed    closed  clamped
remora.ice.bc.ice_thickness      = closed  clamped   closed  clamped
remora.ice.bc.snow_thickness     = closed  closed    closed  closed
remora.ice.bc.meltpond_thickness = closed  closed    closed  closed
remora.ice.bc.ice_age            = closed  closed    closed  closed
remora.ice.bc.Tice               = closed  clamped   closed  clamped
remora.ice.bc.Uice               = closed  clamped   closed  clamped
remora.ice.bc.Vice               = closed  radiation closed  clamped
remora.ice.bc_nudging            = 0 1 0 0

# Zero, so that a variable that does not read the file is seen.
remora.ice.bc_value.Aice          = 0.0 0.0 0.0 0.0
remora.ice.bc_value.ice_thickness = 0.0 0.0 0.0 0.0
remora.ice.bc_value.Tice          = 0.0 0.0 0.0 0.0

# A column away from the corners, which hold the mean of two neighbours.
remora.ice.probe_i = 20
