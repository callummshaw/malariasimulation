#' @title Define model processes
#' @description
#' create_processes, defines the functions which describe how each individual's
#' variables change over time.
#'
#' @param renderer a renderer object
#' @param variables a list of variables in the model
#' @param events a list of events in the model
#' @param parameters a list of model parameters
#' @param models a list of vector models, one for each species
#' @param solvers a list of ode solvers, one for each species
#' @param correlations the intervention correlations object
#' @param lagged_eir a list of list of LaggedValue objects for EIR for each
#' population and species in the simulation
#' @param lagged_infectivity a list of LaggedValue objects for FOIM for each population
#' in the simulation
#' @param timesteps Number of timesteps
#' @param mixing a vector of mixing coefficients for the lagged transmission
#' values (default: 1)
#' @param mixing_fn function to retrieve mixed EIR and infectivity values
#' @param mixing_index an index for this population's position in the
#' lagged transmission lists (default: 1)
#' @noRd
create_processes <- function(
    renderer,
    variables,
    events,
    parameters,
    models,
    solvers,
    correlations,
    lagged_eir,
    lagged_infectivity,
    timesteps,
    mixing_fn = NULL,
    mixing_index = 1
) {
  
  # ========
  # Immunity
  # ========
  processes <- list(
    # Immunity to clinical disease (maternal)
    immunity_process = create_exponential_decay_process(variables$icm,
                                                        parameters$rm),
    # Immunity to clinical disease (acquired)
    immunity_process = create_exponential_decay_process(variables$ica,
                                                        parameters$rc)
  )
  
  if(parameters$parasite == "falciparum"){
    processes <- c(
      processes,
      # Blood immunity
      immunity_process = create_exponential_decay_process(variables$ib,
                                                          parameters$rb),
      # Immunity to severe disease (maternal)
      immunity_process = create_exponential_decay_process(variables$ivm,
                                                          parameters$rvm),
      # Immunity to severe disease (acquired)
      immunity_process = create_exponential_decay_process(variables$iva,
                                                          parameters$rva),
      # Immunity to detectability
      immunity_process = create_exponential_decay_process(variables$id, 
                                                          parameters$rid)
    )
  } else if (parameters$parasite == "vivax"){
    processes <- c(
      processes,
      # Anti-parasite immunity
      immunity_process = create_exponential_decay_process(variables$iam, 
                                                          parameters$rm),
      immunity_process = create_exponential_decay_process(variables$iaa, 
                                                          parameters$ra),
      hypnozoite_process = create_hypnozoite_batch_decay_process(variables$hypnozoites,
                                                                 parameters$gammal)
    )
  }

  if (parameters$individual_mosquitoes) {
    processes <- c(
      processes,
      mosquito_emergence_process = create_mosquito_emergence_process(
        solvers,
        variables$mosquito_state,
        variables$species,
        parameters$species,
        parameters$dpl
      )
    )
  }
  
  # =====================================================
  # Competing Hazard Outcomes (infections and disease progression)
  # =====================================================
  
  if(parameters$parasite == "falciparum"){
    infection_outcome <- CompetingOutcome$new(
      targeted_process = function(timestep, target){
        falciparum_infection_outcome_process(timestep, target, 
                                  variables, renderer, parameters
        )
      },
      size = parameters$human_population
    )
    
  } else if (parameters$parasite == "vivax"){
    infection_outcome <- CompetingOutcome$new(
      targeted_process = function(timestep, target, relative_rates){
        vivax_infection_outcome_process(timestep, target, 
                                  variables, renderer, parameters, 
                                  relative_rates
        )
      },
      size = parameters$human_population
    )
  }
  
  progression_outcome <- CompetingOutcome$new(
    targeted_process = function(timestep, target){
      progression_outcome_process(timestep, target, variables, parameters, renderer)
    },
    size = parameters$human_population
  )
  
  # ==============================
  # Biting and mortality processes
  # ==============================
  # simulate bites, calculates infection rates for bitten humans and set last_boosted_*
  # move mosquitoes into incubating state
  # kill mosquitoes caught in vector control
  processes <- c(
    processes,
    biting_process = create_biting_process(
      renderer,
      solvers,
      models,
      variables,
      events,
      parameters,
      lagged_infectivity,
      lagged_eir,
      mixing_fn,
      mixing_index,
      infection_outcome
    )
  )
  
  # ===================
  # Disease Progression
  # ===================
  
  processes <- c(
    processes,
    progression_process = create_progression_rates_process(
      parameters,
      variables,
      progression_outcome
    ),
    
    # Resolve competing hazards of infection with disease progression
    hazard_resolution_process = CompetingHazard$new(
      outcomes = list(infection_outcome, progression_outcome),
      size = parameters$human_population
    )$resolve
  )
  
  # ===============
  # ODE integration
  # ===============
  processes <- c(
    processes,
    solver_process = create_solver_stepping_process(solvers, parameters)
  )

  # =========
  # PEV EPI
  # =========
  if (!is.null(parameters$pev_epi_coverage)) {
    processes <- c(
      processes,
      epi_pev_process = create_epi_pev_process(
        variables,
        events,
        parameters,
        correlations,
        parameters$pev_epi_coverage,
        parameters$pev_epi_timesteps
      )
    )
  }

  # =========
  # PMC
  # =========
  if(!is.null(parameters$pmc_coverages)){
    processes <- c(
      processes,
      pmc_process = create_pmc_process(
        variables,
        events,
        parameters,
        renderer,
        correlations,
        parameters$pmc_coverages,
        parameters$pmc_timesteps,
        parameters$pmc_drug
      )
    )
  }

  # =========
  # Rendering
  # =========
  
  imm_var_names <- c('ica', 'icm')
  if(parameters$parasite == "falciparum"){
    imm_var_names <- c(imm_var_names, 'ib', 'iva', 'ivm', 'id')
  } else if (parameters$parasite == "vivax"){
    imm_var_names <- c(imm_var_names, 'iaa', 'iam', 'hypnozoites')

    ## hypnozoite infection prevalence rendering
    processes <- c(
      processes,
      hypnozoite_prevalence_renderer = create_n_with_hypnozoites_renderer_process(
        renderer,
        variables$hypnozoites,
        parameters
      ),
      hypnozoite_prevalence_age_group_renderer = create_n_with_hypnozoites_age_renderer_process(
        variables$hypnozoites,
        variables$birth,
        parameters,
        renderer
      )
    )
  }

  processes <- c(
    processes,
    categorical_renderer = individual::categorical_count_renderer_process(
        renderer,
        variables$state,
        c('S', 'A', 'D', 'U', 'Tr')
      ),
    immunity_renderer = create_variable_mean_renderer_process(
      renderer,
      imm_var_names,
      variables[imm_var_names]
    ),
    prevalence_renderer = create_prevalence_renderer(
      variables$state,
      variables$birth,
      variables$id,
      parameters,
      renderer
    ),
    age_group_renderer = create_age_group_renderer(
      variables$birth,
      parameters,
      renderer
    ),
    age_aggregated_immunity_renderer = create_age_variable_mean_renderer_process(
      imm_var_names[paste0(imm_var_names,"_rendering_min_ages") %in% names(parameters)],
      variables[imm_var_names[paste0(imm_var_names,"_rendering_min_ages") %in% names(parameters)]],
      variables$birth,
      parameters,
      renderer
    ),
    mosquito_state_renderer = create_compartmental_rendering_process(
      renderer,
      solvers,
      parameters
    )
  )

  if (parameters$individual_mosquitoes) {
    processes <- c(
      processes,
      vector_count_renderer = create_vector_count_renderer_individual(
        variables$mosquito_state,
        variables$species,
        variables$mosquito_state,
        renderer,
        parameters
      )
    )
  } else {
    processes <- c(
      processes,
      vector_count_renderer = create_total_M_renderer_compartmental(
        renderer,
        solvers,
        parameters
      )
    )
  }

  # ======================
  # Intervention processes
  # ======================

  if (parameters$bednets) {
    processes <- c(
      processes,
      distribute_nets_process = distribute_nets(
        variables,
        events$throw_away_net,
        events$change_net,
        parameters,
        correlations
      ),
      net_usage_renderer = net_usage_renderer(variables$net_time, renderer)
    )
  }

  if (parameters$spraying) {
    processes <- c(
      processes,
      indoor_spraying_process = indoor_spraying(
        variables$spray_time,
        renderer,
        parameters,
        correlations
      )
    )
  }

  # ======================
  # Progress bar process
  # ======================
  if (parameters$progress_bar){
    processes <- c(
      processes,
      progress_bar_process = create_progress_process(timesteps)
    )
  }

  # ======================
  # Mortality step
  # ======================
  # Mortality is not resolved as a competing hazard
  
  processes <- c(
    processes,
    mortality_process = create_mortality_process(
      variables,
      events,
      renderer,
      parameters
    )
  )

  processes
}

# =================
# Utility functions
# =================

#' @title Exponentially decaying variables
#' @description
#' create_exponential_decay_process generates a process function
#' that reduces the value of a variable at an exponential rate
#'
#' @param variable the variable to update
#' @param rate the exponential rate
#' @noRd
create_exponential_decay_process <- function(variable, rate) {
  stopifnot(inherits(variable, "DoubleVariable"))
  decay_rate <- exp(-1/rate)
  exponential_process_cpp(variable$.variable, decay_rate)
}

#' @title Create and initialise lagged_infectivity object
#'
#' @param variables model variables for initialisation
#' @param parameters model parameters
#' @noRd
create_lagged_infectivity <- function(variables, parameters) {
  age <- get_age(variables$birth$get_values(), 0)
  psi <- unique_biting_rate(age, parameters)
  .pi <- human_pi(psi, variables$zeta$get_values())
  init_infectivity <- sum(.pi * variables$infectivity$get_values())
  LaggedValue$new(
    max_lag = parameters$delay_gam + 2,
    default = init_infectivity
  )
}

#' @title Create and initialise a list of lagged_eir objects per species
#'
#' @param variables model variables for initialisation
#' @param solvers model differential equation solvers
#' @param parameters model parameters
#' @noRd
create_lagged_eir <- function(variables, solvers, parameters) {
  lapply(
    seq_along(parameters$species),
    function(species) {
      LaggedValue$new(
        max_lag = parameters$de + 1,
        default = calculate_eir(
          species,
          solvers,
          variables,
          parameters,
          0
        )
      )
    }
  )
}
#' @title Hypnozoite decay function
#' @description
#' calulates the number of individuals in whom a batch decay occurs
#'
#' @param variable the hypnozoite variable to update
#' @param rate the hypnozoite decay rate
#' @noRd
create_hypnozoite_batch_decay_process <- function(hypnozoites, gammal){
  function(timestep){
    to_decay <- bernoulli_multi_p(p = rate_to_prob(hypnozoites$get_values() * gammal))
    if(length(to_decay) > 0){
      hypnozoites$queue_update(
        hypnozoites$get_values(to_decay) - 1,
        to_decay
      )
    }
  }
}