
# SeaIce_stress_radiation: the internal ice stresses radiate through the west side of the
# thick-ice case, where every ice boundary is otherwise closed.
#                       WEST      SOUTH    EAST   NORTH
remora.ice.bc.ice_Sxx = radiation periodic closed periodic
remora.ice.bc.ice_Sxy = radiation periodic closed periodic
remora.ice.bc.ice_Syy = radiation periodic closed periodic

# The probe is the cell on the radiating side, whose along-shore momentum reads the stress
# the condition sets. It is read after three steps.
remora.ice.probe_i = 0
remora.ice.probe_j = 2
remora.ice.probe_int = 1
max_step = 3
