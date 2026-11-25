
# SeaIce_layout_ns: SeaIce_BC with the south and north sides open instead of the west and east.
# The layout test cuts these sides along their length, so that the ice in the row past them is
# read across a box edge.
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
remora.ice.bc.Aice               = closed  clamped   closed  radiation
remora.ice.bc.ice_thickness      = closed  clamped   closed  radiation
remora.ice.bc.snow_thickness     = closed  closed    closed  radiation
remora.ice.bc.meltpond_thickness = closed  closed    closed  radiation
remora.ice.bc.ice_age            = closed  closed    closed  radiation
remora.ice.bc.ice_Sxx            = closed  closed    closed  radiation
remora.ice.bc.ice_Sxy            = closed  closed    closed  radiation
remora.ice.bc.ice_Syy            = closed  closed    closed  radiation
remora.ice.bc.Tice               = closed  clamped   closed  gradient
remora.ice.bc.Uice               = closed  clamped   closed  radiation
remora.ice.bc.Vice               = closed  clamped   closed  radiation
remora.ice.bc_nudging            = 0 0 0 1

remora.ice.bc_value.Aice          = 0.0 0.5  0.0 0.0
remora.ice.bc_value.ice_thickness = 0.0 0.2  0.0 0.0
remora.ice.bc_value.Tice          = 0.0 -4.0 0.0 0.0
remora.ice.bc_value.Uice          = 0.0 0.1  0.0 0.0
remora.ice.bc_value.Vice          = 0.0 0.1  0.0 0.0
