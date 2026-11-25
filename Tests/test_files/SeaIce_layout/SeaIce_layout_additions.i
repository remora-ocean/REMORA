
# SeaIce_layout: SeaIce_BC on a grid wide enough to cut its periodic direction, with walls for
# the ocean so that only the ice has open sides.
amr.n_cell = 39 16 30
geometry.prob_hi = 100000. 20000. 0.
remora.bc.temp.type   = slipwall periodic slipwall periodic
remora.bc.salt.type   = slipwall periodic slipwall periodic
remora.bc.u.type      = slipwall periodic slipwall periodic
remora.bc.v.type      = slipwall periodic slipwall periodic
remora.bc.ubar.type   = slipwall periodic slipwall periodic
remora.bc.vbar.type   = slipwall periodic slipwall periodic
remora.bc.zeta.type   = slipwall periodic slipwall periodic
remora.bc.tke.type    = slipwall periodic slipwall periodic
