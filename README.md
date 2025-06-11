[![R build status](https://github.com/mrc-ide/malariasimulation/workflows/R-CMD-check/badge.svg)](https://github.com/mrc-ide/malariasimulation/actions)
[![codecov](https://codecov.io/github/mrc-ide/malariasimulation/branch/master/graphs/badge.svg)](https://codecov.io/github/mrc-ide/malariasimulation)

# malariasimulation <img src="man/figures/malariasimulation.png" align="right" width=30% height=30% />

A Fork of Imperial College London's next malaria simulation. I have added multiple ITN functionality with seperate decays for rates for repellency and mortality. 

Nets are distributed at the same timestep but the proportion of each net type distributed can be varied in time, similarily other attributes can vary time, net type, and vector species 
The bednets are now set like --

Assuming:
nt - number of net distribution times
nv - number of vector species
nn - number of net type 

model_input <- set_bednets(
    sim_params,
    timesteps = times, #Times we want to distribute the nets (nt x 1) 
    coverages = net_params$cov, #The overall coverages for each distribution time, proportion of pop that gets ANY net (nt x 1)
    coverage_by_type = net_params$cov_n, #The relative distribution of each net type, can vary with time (nt x nn)
    retention = net_params$ret, #How long nets are kept for currenly a constant but would be able to make vary by time or net type 
    dn0 = net_params$dn0, #Initial death from nets (nt x nv x nn)
    rn = net_params$rn,  #Initial repellency from nets (nt x nv x nn)
    rnm = net_params$rnm,  #Base repellency from nets (nt x nv x nn)
    gammad =  net_params$gammad, #Death decay rate (nt x nn)
    gammar =  net_params$gammar #Repellency decay rate (nt x nn)
  )

