

population_grouping <- structure(list(`Young Male Black` = 1:2, `Old Male Black` = 3:4,
                                      `Young Female Black` = 17:18, `Old Female Black` = 19:20,
                                      `Young Male Hispanic` = 9:10, `Old Male Hispanic` = 11:12,
                                      `Young Female Hispanic` = 25:26, `Old Female Hispanic` = 27:28,
                                      `Young Male Other` = 5:6, `Old Male Other` = 7:8, `Young Female Other` = 21:22,
                                      `Old Female Other` = 23:24, `Young Male MSM` = 13:14, `Old Male MSM` = 15:16,
                                      `Young Female MSM` = 29:30, `Old Female MSM` = 31:32))


# 1+d.index are the diagnosed columns in the sol data frame.
d.index <- s.index + 32*4

# Incidence index
i.index <- s.index + 32*6

# Screening Index
scr_index <- s.index + 32*7

# Relevant Data Rows
# change from 2000:2015 -> 2010:2021, i.e. add 10 and 6
model_years <- 2010:2021
model_year_rows <- 56:67
model_year_rows_1 <- 55:67

get_screening_rate_from_population <- function(popname, theta, abcd_point = 'a') {
  sex <-
    gsub("Black|Hispanic|Other|MSM|Young|Old| ", "", popname)
  demographic <-
    gsub("Young|Old|Female|Male| ", "", popname)
  age <-
    gsub("Black|Hispanic|Other|MSM|Female|Male| ", "", popname)

  # attach abcd_point to the end of string
  p <- function(a) paste0(a, abcd_point)

  if (demographic %in% c('Hispanic', 'Other') && sex == 'Male' && age == "Young") param <- p('logit.screen.m1.')
  if (demographic %in% c('Hispanic', 'Other') && sex == 'Male' && age == "Old") param <- p('logit.screen.m2.')
  if (demographic=='Black' && sex == 'Male' && age == "Young") param <- p('logit.screen.m1.1.')
  if (demographic=='Black' && sex == 'Male' && age == "Old") param <- p('logit.screen.m2.1.')
  if (demographic=='MSM' && age == "Young") param <- p('logit.screen.msm1.')
  if (demographic=='MSM' && age == "Old") param <- p('logit.screen.msm2.')
  if (demographic %in% c('Hispanic', 'Other') && sex == 'Female' && age == "Young") param <- p('logit.screen.f1.')
  if (demographic %in% c('Hispanic', 'Other') && sex == 'Female' && age == "Old") param <- p('logit.screen.f2.')
  if (demographic=='Black' && sex == 'Female' && age == "Young") param <- p('logit.screen.f1.1.')
  if (demographic=='Black' && sex == 'Female' && age == "Old") param <- p('logit.screen.f2.1.')

  param <- theta[[param]]
  param <- ilogit(param)

  # I think we need to multiply by RR screen for Hispanic, and somehow account
  # for the RR screening of High Activity

  return(param)
}

# Get a dataframe out of sol which has columns:
# year, population, population_size, symptomatic_cases, asymptomatic_cases, diagnosed_cases
extract_data_from_sol <- function(sol) {
  simulation_data <- list()
  # Loop over Population Groups
  for (i in seq(population_grouping)) {
    # Get their population indices
    pop <- population_grouping[[i]]
    # Get the population name
    popname <- names(population_grouping)[[i]]
    # Number symptomatic
    symptomatic_data <- rowSums(sol[model_year_rows, 1+y.index[pop]])
    # Number asymptomatic
    asymptomatic_data <- rowSums(sol[model_year_rows, 1+z.index[pop]])
    # Total population size
    totalpop_data <- symptomatic_data + asymptomatic_data + rowSums(sol[model_year_rows, 1+s.index[pop]])
    # Use cumulative diagnoses and subtract the next year from the previous to get yearly diagnoses
    diagnosed_data <- rowSums(sol[model_year_rows_1, 1+d.index[pop]])
    diagnosed_data <- diagnosed_data[2:length(diagnosed_data)] - diagnosed_data[1:(length(diagnosed_data)-1)]
    # Use cumulative incidence and subtract the next year from the previous to get yearly incidence
    incident_data <- rowSums(sol[model_year_rows_1, 1+i.index[pop]])
    incident_data <- incident_data[2:length(incident_data)] - incident_data[1:(length(incident_data)-1)]
    # Use cumulative screened cases and subtract the next year from the previous to get yearly screened cases
    screened_data <- rowSums(sol[model_year_rows_1, 1+scr_index[pop]])
    screened_data <- screened_data[2:length(screened_data)] - screened_data[1:(length(screened_data)-1)]

    # Construct population specific dataframe
    simulation_data[[length(simulation_data)+1]] <-
      cbind.data.frame(
        years = model_years,
        population = popname,
        population_size = totalpop_data,
        symptomatic_cases = symptomatic_data,
        asymptomatic_cases = asymptomatic_data,
        diagnosed_cases = diagnosed_data,
        incident_cases = incident_data,
        number_screened = screened_data
      )
  }
  # Rbind all dataframes together.
  return(do.call(rbind.data.frame, simulation_data))
}


# Trim down the trace.
# trace <- trace[100:nrow(trace), ]
# trace <- dplyr::sample_n(trace, 100)

duration_param_names <-
  c("log.dur.inf.symp.m",
    "log.dur.inf.symp.msm",
    "log.dur.inf.symp.f",
    "log.dur.inf.asymp.m",
    "log.dur.inf.asymp.msm",
    "log.dur.inf.asymp.f")

pr_symp_param_names <-
  c("logit.symp.m",
    "logit.symp.msm",
    "logit.symp.f")


insert_duration_and_symptomatic_param_cols <- function(theta, sol) {

  duration_params <- theta[duration_param_names]
  pr_symp_params <- theta[pr_symp_param_names]

  duration_params <- exp(duration_params)
  pr_symp_params <- ilogit(pr_symp_params)

  names(duration_params) <- gsub("log.", "", names(duration_params))
  names(pr_symp_params) <- gsub("logit.", "", names(pr_symp_params))

  # Get duration of symptomatic cases
  duration_symptomatic_param <- sapply(sol$population, function(pop) {
    if (grepl("Male.*MSM", pop)) { # msm
      return(duration_params[['dur.inf.symp.msm']])
    } else if (grepl("Male", pop)) { # heterosexual male
      return(duration_params[['dur.inf.symp.m']])
    } else { # female
      return(duration_params[['dur.inf.symp.f']])
    }
  })

  # Get duration of symptomatic cases
  duration_asymptomatic_param <- sapply(sol$population, function(pop) {
    if (grepl("Male.*MSM", pop)) { # msm
      return(duration_params[['dur.inf.asymp.msm']])
    } else if (grepl("Male", pop)) { # heterosexual male
      return(duration_params[['dur.inf.asymp.m']])
    } else { # female
      return(duration_params[['dur.inf.asymp.f']])
    }
  })

  # Get probability symptomatic
  pr_symptomatic_case_param <- sapply(sol$population, function(pop) {
    if (grepl("Male.*MSM", pop)) { # msm
      return(pr_symp_params[['symp.msm']])
    } else if (grepl("Male", pop)) { # heterosexual male
      return(pr_symp_params[['symp.m']])
    } else { # female
      return(pr_symp_params[['symp.f']])
    }
  })

  # Insert duration and symptomatic probabilities into sol
  sol <- cbind.data.frame(
    sol,
    prob_symptomatic_case_param = pr_symptomatic_case_param,
    duration_symptomatic_param = duration_symptomatic_param,
    duration_asymptomatic_param = duration_asymptomatic_param
  )

  return(sol)
}



simulate_and_save_simulation_stats_and_ci_stats <- function(filepath){
  df <- foreach (iter=1:nrow(trace), .combine = rbind.data.frame) %do% {
    try({
      # get theta and run simulation
      theta <- unlist(trace[iter,])
      e <- create_gcSim_param_env(theta, interv = interv)
      gc_env$e <- e
      run_gcSim_with_environment(e)
      sol <- e$out.cpp$out
      sol <- cbind.data.frame(interv = interv, simulation.n = iter, extract_data_from_sol(sol))
      sol <- insert_duration_and_symptomatic_param_cols(theta, sol)
    })
  }

  # add screening rates ----
  screen_a <- numeric(length = nrow(df))
  foreach (i=seq(nrow(df))) %do% {
    theta <- unlist(trace[df[[i,'simulation.n']], ])

    screen_a[[i]] <- get_screening_rate_from_population(
      popname = df[[i, 'population']],
      theta = theta,
      abcd_point = 'a'
    )
  }
  df$screen_a <- screen_a

  screen_d <- numeric(length = nrow(df))
  foreach (i=seq(nrow(df))) %do% {
    theta <- unlist(trace[df[[i,'simulation.n']], ])

    screen_d[[i]] <- get_screening_rate_from_population(
      popname = df[[i, 'population']],
      theta = theta,
      abcd_point = 'd'
    )
  }
  df$screen_d <- screen_d

  rr.screen.hispanic <- numeric(length = nrow(df))
  for (i in seq(nrow(df))) {
    theta <- unlist(trace[df[[i,'simulation.n']], ])
    if (grepl("Male Hispanic", df[[i, 'population']])) {
      rr.screen.hispanic[[i]] <- exp(theta[['log.rr.screen.m3']])
    } else if (grepl("Female Hispanic", df[[i, 'population']])) {
      rr.screen.hispanic[[i]] <- exp(theta[['log.rr.screen.f3']])
    } else {
      rr.screen.hispanic[[i]] <- 1
    }
  }
  df$rr.screen.hispanic <- rr.screen.hispanic

  rr.high_ac.screen <- numeric(length = nrow(df))
  for (i in seq(nrow(df))) {
    theta <- unlist(trace[df[[i,'simulation.n']], ])
    rr.high_ac.screen[[i]] <- exp(theta[['log.rr.screen.ac']])
  }
  df$rr.screen.high_ac <- rr.high_ac.screen


  # clean data ----
  # Get rid of Female MSM
  df <- filter(df, ! grepl(".* Female MSM", population))
  # Add prevalence, proportion symptomatic, and proportion of incidence diagnosed.
  df <- mutate(df,
               prevalence = (symptomatic_cases + asymptomatic_cases) / population_size,
               proportion_symptomatic = symptomatic_cases / (symptomatic_cases + asymptomatic_cases),
               proportion_diagnosed = diagnosed_cases / incident_cases,
               proportion_tested = number_screened / population_size)
  # Let's filter out the simulations where any negative numbers of symptomatic cases appear.
  bad_theta_idxs <- filter(df, symptomatic_cases < 0 | number_screened < 0) %>% `[[`('simulation.n')
  df <- filter(df, ! simulation.n %in% bad_theta_idxs)
  # Rename number_screened -> number_tested and
  colnames(df)[[which(colnames(df)=='number_screened')]] <- 'number_tested'

  # save first dataframe ----
  saveRDS(object = df, file = paste0(filepath, interv, "/simulation statistics.rds"))

  # summarize and save second data frame ----
  dt <-
    select(
      df,
      years,
      population,
      population_size,
      prevalence,
      proportion_symptomatic,
      proportion_diagnosed,
      proportion_tested
    )

  # reshape for plotting
  dt <- reshape2::melt(dt, id.vars = c('years', 'population'))

  # summarize by population in each year for each statistical variable
  dt <- dt %>% group_by(years, population, variable) %>%
    summarise(mean = mean(value, na.rm = T),
              ci_low = quantile(value, .025, na.rm = T),
              ci_high = quantile(value, .975, na.rm = T))

  # arrange and save
  dt <- arrange(dt, population, variable, years)

  write.csv(dt, file = paste0(filepath, interv, "/confidence intervals stats.csv"))
}




#' Prevalence Estimates.
#'
#' We want to output Prevalence by
#'     - Gender, race-ethnicity, (grouping risk, age)
#'     - Race-ethnicity, (grouping gender, risk, age)
#'     - Overall prevalence
#'     - Gender (grouping race-ethnicity, risk, age)
#'     - Gender, age (grouping race-ethnicity, risk)
#'
#' @output A dataframe with the following columns:
#'     - year
#'     - total_population_size
#'     - overall
#'     - black_male
#'     - hispanic_male
#'     - other_male
#'     - black_female
#'     - hispanic_female
#'     - other_female
#'     - msm
#'     - young msm
#'     - old msm
#'     - black
#'     - hispanic
#'     - other
#'     - male
#'     - female
#'     - male_young
#'     - male_old
#'     - female_young
#'     - female_old
#' All columns except `year` and `total_population_size` will be prevalence rates.
#'
#' @example
#' source("Model_code/load.start.conditions.R")
#' load("Model_outputs/Main/cal/merged_trace")
#' source("~/Documents/2018/05 May 2018/May29/gcNational_prevalence_outputs_for_disparities.R")
#' load.start()
#' theta <- trace[sample.int(n = nrow(trace), size=1), ]
#' e <- construct_prediction_environment(theta)
#' out.cpp <- run_simulation_with_env(e)
#' sol <- out.cpp$out
#' prevalence_estimates_for_disparities(sol)
prevalence_estimates_for_disparities <- function(sol) {
  df_colnames <-
    c(
      'year',
      'overall',
      'black_male',
      'hispanic_male',
      'other_male',
      'black_female',
      'hispanic_female',
      'other_female',
      'msm',
      'young msm',
      'old msm',
      'black',
      'hispanic',
      'other',
      'male',
      'female',
      'male_young',
      'male_old',
      'female_young',
      'female_old'
    )

  # vector of years which we will measure
  year <- 2010:2021
  year_idx <- (cal.start+10):(81)

  df <- list()
  df$year <- year

  all_infections <-
    sol[year_idx,1+y.index] +
    sol[year_idx,1+z.index]

  population_sizes <-
    all_infections +
    sol[year_idx,1+s.index] +
    sol[year_idx,1+nsa.index]

  # get rid of the "Female MSM"
  all_infections <- all_infections[,1:28]
  population_sizes <- population_sizes[,1:28]

  # calculating overall prevalence
  overall <- rowSums(all_infections)/rowSums(population_sizes)
  df$overall <- overall

  population_prevalence <- function(idx) rowSums(all_infections[, idx]) / rowSums(population_sizes[, idx])

  # save prevalence rates for various groups
  df$black_male <- population_prevalence(m1)
  df$hispanic_male <- population_prevalence(m3)
  df$other_male <- population_prevalence(m2)
  df$msm <- population_prevalence(m4)
  df$young_msm <- population_prevalence(y.m.4)
  df$old_msm <- population_prevalence(o.m.4)
  df$black_female <- population_prevalence(f1)
  df$hispanic_female <- population_prevalence(f3)
  df$other_female <- population_prevalence(f2)
  df$black <- population_prevalence(c(m1,f1))
  df$hispanic <- population_prevalence(c(m3,f3))
  df$other <- population_prevalence(c(m2,f2))
  df$male <- population_prevalence(c(m1,m2,m3,m4))
  df$female <- population_prevalence(c(f1,f2,f3))
  df$male_young <- population_prevalence(y.m)
  df$male_old <- population_prevalence(o.m)
  df$female_young <- population_prevalence(y.f)
  df$female_old <- population_prevalence(o.f)

  # return output dataframe
  df <- cbind.data.frame(df)
  return(df)
}


#' Prevalence Estimates with Confidence Intervals for Disparities
#'
#' This function will take a dataframe a sample of a trace, run each
#' parameter vector as a simulation, use prevalence_estimates_for_disparities
#' on the simulation data, melt the data, and store it in a list of data
#' frames which are then used to compute confidence intervals.
#' @example
#' df <- prevalence_estimates_with_confidence_for_disparities(trace[1:3,])
#' library(ggplot2)
#' ggplot(filter(df, group!='msm'), # plot prevalence for all groups except MSM
#'        aes(x=year,
#'            y=mean,
#'            ymax = ci_high,
#'            ymin = ci_low,
#'            group = group,
#'            color = group,
#'            fill = group)) +
#'   geom_line() +
#'   geom_ribbon(alpha=0.4, size=0)
prevalence_estimates <- function(trace) {

  extract_theta <- function(i) {
    theta <- as.numeric(trace[i,])
    names(theta) <- colnames(trace)
    return(theta)
  }

  all_sims_prev <- list()
  for(iter in 1:nrow(trace)) {
    theta <- extract_theta(iter)
    e <- create_gcSim_param_env(theta, interv)
    run_gcSim_with_environment(e)
    df <- prevalence_estimates_for_disparities(e$out.cpp$out)
    df <- reshape2:::melt.data.frame(df, id.vars = c('year'), variable.name = 'group')
    df <- cbind.data.frame(sim = iter, df)
    all_sims_prev[[iter]] <- df
  }
  # all_sims_prev <- rbind(all_sims_prev)
  all_sims_prev <- do.call(rbind, all_sims_prev)

  return(all_sims_prev)
}
