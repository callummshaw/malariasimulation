#' @title Parameterise a bed net strategy
#'
#' @description The model will distribute bed nets at `timesteps` to a random
#' sample of the entire human population. The sample size will be a proportion
#' of the human population taken from the corresponding `coverages`.
#' The sample _can_ contain humans who already have bed nets.
#'
#' All of the sample "use" their bed nets on the timestep after they are
#' distributed. Incomplete usage is not part of this model.
#'
#' If a human in the sample already has a bed net, their bed net will be replaced
#' by a new one.
#'
#' Bed nets will be randomly removed each timestep with a rate of `1 -
#' exp(-1/retention)`
#'
#' The structure for the bed net model is documented in the
#' S.I. of 10.1038/s41467-018-07357-w
#'
#' @param parameters a list of parameters to modify
#' @param timesteps the timesteps at which to distribute bed nets
#' @param coverages the proportion of the human population who receive bed nets
#' @param retention the average number of timesteps a net is kept for (an array of the number of nets)
#' @param dn0 a matrix of death probabilities for each species over time.
#' With nrows=length(timesteps), ncols=length(species)
#' @param rn a matrix of repelling probabilities for each species over time
#' With nrows=length(timesteps), ncols=length(species)
#' @param rnm a matrix of minimum repelling probabilities for each species over time
#' With nrows=length(timesteps), ncols=length(species)
#' @param gammad a vector of bednet killing half-lives for each distribution timestep
#' @param gammar a vector of bednet repellency half-lives for each distribution timestep

#' @export
set_bednets <- function(
    parameters,
    timesteps,
    coverages,
    coverage_by_type,
    retention,
    dn0,
    rn,
    rnm,
    gammad,
    gammar,
    replace
) {
  stopifnot(all(coverages >= 0) && all(coverages <= 1))
  
  parameters$bednets <- TRUE
  parameters$bednet_timesteps <- timesteps
  parameters$bednet_coverages <- coverages
  parameters$bednet_coverage_by_type <- coverage_by_type
  parameters$bednet_dn0 <- dn0
  parameters$bednet_rn <- rn
  parameters$bednet_rnm <- rnm
  parameters$bednet_gammad <- gammad
  parameters$bednet_gammar <- gammar
  parameters$bednet_retention <- retention
  parameters$bednet_replace <- replace
  parameters
}

#' @title Parameterise an indoor spraying strategy
#'
#' @description The model will apply indoor spraying at `timesteps` to a random
#' sample of the entire human population. The sample size will be a proportion
#' of the human population taken from the corresponding `coverages`.
#' The sample _can_ contain humans who have already benefited from spraying.
#'
#' If a human in the sample lives in a sprayed house, the efficacy of the
#' spraying will be returned to the maximum.
#'
#' The structure for the indoor residual spraying model is documented in the
#' S.I. of 10.1038/s41467-018-07357-w
#'
#' @param parameters a list of parameters to modify
#' @param timesteps the timesteps at which to spray
#' @param coverages the proportion of the population who get indoor
#' spraying
#' @param ls_theta matrix of mortality parameters
#' With nrows=length(timesteps), ncols=length(species)
#' @param ls_gamma matrix of mortality parameters per timestep
#' With nrows=length(timesteps), ncols=length(species)
#' @param ks_theta matrix of feeding success parameters per timestep 
#' With nrows=length(timesteps), ncols=length(species)
#' @param ks_gamma matrix of feeding success parameters per timestep 
#' With nrows=length(timesteps), ncols=length(species)
#' @param ms_theta matrix of deterrence parameters per timestep 
#' With nrows=length(timesteps), ncols=length(species)
#' @param ms_gamma matrix of deterrence parameters per timestep 
#' With nrows=length(timesteps), ncols=length(species)
#' @export
set_spraying <- function(
    parameters,
    timesteps,
    coverages,
    ls_theta,
    ls_gamma,
    ks_theta,
    ks_gamma,
    ms_theta,
    ms_gamma
) {
  stopifnot(all(coverages >= 0) && all(coverages <= 1))
  if (length(coverages) != length(timesteps)) {
    stop('coverages and timesteps must must align')
  }
  decays <- list(
    ls_theta,
    ls_gamma,
    ks_theta,
    ks_gamma,
    ms_theta,
    ms_gamma
  )
  for (x in decays) {
    if (ncol(x) != length(parameters$species)) {
      stop('theta and gamma rows need to align with species')
    }
    if (nrow(x) != length(timesteps)) {
      stop('theta and gamma cols need to align with timesteps')
    }
  }
  parameters$spraying <- TRUE
  parameters$spraying_timesteps <- timesteps
  parameters$spraying_coverages <- coverages
  parameters$spraying_ls_theta <- ls_theta
  parameters$spraying_ls_gamma <- ls_gamma
  parameters$spraying_ks_theta <- ks_theta
  parameters$spraying_ks_gamma <- ks_gamma
  parameters$spraying_ms_theta <- ms_theta
  parameters$spraying_ms_gamma <- ms_gamma
  parameters
}

#' @title Parameterise custom baseline carrying capacity
#' 
#' @description Allows the user to set a completely flexible and custom
#' carrying capacity for each species
#' 
#' @param parameters the model parameters
#' @param timesteps vector of timesteps for each rescale change
#' @param carrying_capacity_scalers matrix of scaling factors to scale the baseline 
#' carrying capacity for each species with nrows = length(timesteps),
#'  ncols = length(species)
#' 
#' @export
set_carrying_capacity <- function(
    parameters,
    timesteps,
    carrying_capacity_scalers
){
  stopifnot(nrow(carrying_capacity_scalers) == length(timesteps))
  stopifnot(ncol(carrying_capacity_scalers) == length(parameters$species))
  stopifnot(min(timesteps) > 0)
  stopifnot(min(carrying_capacity_scalers) >= 0)
  
  parameters$carrying_capacity <- TRUE
  parameters$carrying_capacity_timesteps <- timesteps
  parameters$carrying_capacity_scalers <- carrying_capacity_scalers
  parameters
}

#' Get initialised carrying capacity for each species
#'
#' @param parameters the model parameters
#'
#' @return a vector of carrying initialised carrying capacity estimates for
#' each vector species 
#' @export
get_init_carrying_capacity <- function(parameters){
  init_cc <- sapply(1:length(parameters$species), function(x){
    p <- parameters$species_proportions[[x]]
    m <- p * parameters$total_M
    calculate_carrying_capacity(parameters, m, x)
  })
  names(init_cc) <- parameters$species
  return(init_cc)
}
