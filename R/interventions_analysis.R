#########################
# Intervention Analysis #
#########################


#' Consruct the List of Key Populations which will be used in the intervention analysis.
#'
#' This function is called in load_start so that the key_populations will be available
#' for the intervention_analysis.
construct_key_populations <- function() {
  stopifnot(exists('gc_env'))
  with(gc_env, {
  key_populations <- list(
    all = 1:32,
    black_male = m1,
    hispanic_male = m3,
    other_male = m2,
    msm = m4,
    msm_young = y.m.4,
    msm_old = o.m.4,
    black_female = f1,
    hispanic_female = f3,
    other_female = f2,
    black = c(m1, f1),
    hispanic = c(m3, f3),
    other = c(m2, f2),
    male = c(m1, m2, m3, m4),
    female = c(f1, f2, f3),
    male_young = y.m,
    male_old = o.m,
    female_young = y.f,
    female_old = o.f
  )

  key_population_names <- list(
    all = "All",
    black_male  = "Black Male",
    hispanic_male  = "Hispanic Male",
    other_male  = "Other Male",
    msm = "MSM",
    msm_young  = "MSM Young",
    msm_old  = "MSM Old",
    black_female  = "Black Female",
    hispanic_female  = "Hispanic Female",
    other_female  = "Other Female",
    black  = "Black",
    hispanic  = "Hispanic",
    other  = "Other",
    male  = "Male",
    female  = "Female",
    male_young = "Male Young",
    male_old  = "Male Old",
    female_young = "Female Young",
    female_old = "Female Old"
  )
  })
}

#' Construct Key Metrics from a Simulation out of gcSim
#'
#' Let's take the simulation data from gcSim and save only some key metrics
#' in a data frame.
#'
#'
intervention_analysis_key_metrics <- function(sol) {
  gc_env$sol <- sol
  with(gc_env, {

  simulation_data <- list()
  years_for_output <- c(min(cal_and_intervention_rows)-1, cal_and_intervention_rows)

  # Filter the solution for the years we want
  rows_we_want <- which(sol[,1] %in% years_for_output)
  rows_we_want <- c(min(rows_we_want)-1, rows_we_want)
  sol <- sol[rows_we_want, ]

  sol[,1] <- sol[,1] + (sol[,2]/13)*.25
  sol <- sol[,-2]

  # Loop over Population Groups
  for (iter in seq(key_populations)) {
    # Get their population indices
    pop <- key_populations[[iter]]
    # Get the population name
    popname <- key_population_names[[iter]]
    # Number symptomatic
    symptomatic_data <- rowSums(sol[, 1+y.index[pop]])
    # Number asymptomatic
    asymptomatic_data <- rowSums(sol[, 1+z.index[pop]])
    # Total population size
    totalpop_data <- symptomatic_data + asymptomatic_data + rowSums(sol[, 1+s.index[pop]])
    # Use cumulative diagnoses and subtract the next year from the previous to get yearly diagnoses
    diagnosed_data <- rowSums(sol[, 1+diag.index[pop]])
    # diagnosed_data <- diagnosed_data[2:length(diagnosed_data)] - diagnosed_data[1:(length(diagnosed_data)-1)]
    # Use cumulative incidence and subtract the next year from the previous to get yearly incidence
    incident_data <- rowSums(sol[, 1+inc.index[pop]])
    # incident_data <- incident_data[2:length(incident_data)] - incident_data[1:(length(incident_data)-1)]
    # Use cumulative screened cases and subtract the next year from the previous to get yearly screened cases
    screened_data <- rowSums(sol[, 1+scr.index[pop]])
    # screened_data <- screened_data[2:length(screened_data)] - screened_data[1:(length(screened_data)-1)]

    # Construct population specific dataframe
    simulation_data[[length(simulation_data)+1]] <-
      cbind.data.frame(
        year = sol[,1]+1942,
        population = popname,
        population_size = totalpop_data,
        symptomatic_cases = symptomatic_data,
        asymptomatic_cases = asymptomatic_data,
        diagnosed_cases = diagnosed_data,
        incident_cases = incident_data,
        number_screened = screened_data
      )
  }})
  # Rbind all dataframes together.
  return(do.call(rbind.data.frame, gc_env$simulation_data))
}


melt_intervention_df <- function(intervention_outcomes) {
  reshape2::melt(intervention_outcomes, id.vars = c('year', 'population'))
}

#' Simulate and Interventions for a Trace
simulate_interventions <- function(interventions, trace) {
	if (missing(interventions) || is.null(interventions)) {
		interventions <- c(
			'base_case',
			'no_screening',
			'universal_annual',
			'universal_semi_annual',
			'high_activity_young_annual',
			'high_activity_young_semi_annual',
			'young_annual',
			'young_semi_annual',
			'high_activity_MSM_annual',
			'high_activity_MSM_semi_annual',
			'high_activity_MSM_quarter_annual',
			"MSM_annual",
			"MSM_semi_annual",
			"MSM_quarter_annual",
			"high_activity_annual",
			"high_activity_semi_annual",
			"female_young_annual",
			"female_young_semi_annual",
			"high_activity_female_young_annual",
			"high_activity_female_young_semi_annual",
			"high_activity_female_young_quarter_annual",
			"remove_ltfu_10pct",
			"remove_ltfu_20pct",
			"remove_f_msm_10pct_and_msw_20pct_ltfu",
			"1year_high_activity_semi_annual"
		)
	}

  if (missing(trace)) { 
    trace <- switch(
      site,
      BA = as.data.frame(readRDS(system.file("calibration_outcomes/BA_posterior_sample.rds", package = "gcRegional"))$theta.list),
      SF = as.data.frame(readRDS(system.file("calibration_outcomes/SF_posterior_sample.rds", package = "gcRegional"))$theta.list)
    )
  }

  key_metrics <- list()

  for (iter in 1:nrow(trace)) {
    theta <- unlist(trace[iter, ])
    # simulate each intervention for a parameter vector
    for (interv in interventions) {
      e <- create_gcSim_param_env(theta, interv = interv)
      e$quarterly <- T
      run_gcSim_with_environment(e)
      sol <- e$out.cpp$out
      df <- intervention_analysis_key_metrics(sol)
      key_metrics[[length(key_metrics)+1]] <- cbind.data.frame(sim = iter, interv = interv, df)
      rm(e)
    }
  }
  results <- do.call(rbind, key_metrics)
  return(results)
}
