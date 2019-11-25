#' Simulate Each Theta from a Trace and Calculate Important Values
#'
#' We want the following measures:
#'
#' Cases Averted
#' Prevalence
#' Number of Tests
#' Number of Tests to Avert One Case
#'
#' I think I'll just make a nested list, containing the
#' cpp outcomes. The top layer will be integer indexed
#' in coordination with a calibration trace. Inside
#' each single parameter vector from the trace will
#' correspond to a named list of simulations.
#'
# simulate_interventions <- function(site, years, partial_results_path) {
#   load_start(site)
#
#   interventions <- c(
#     'base_case',
#     'high_activity_young_annual',
#     'high_activity_young_semi_annual',
#     'young_annual',
#     'young_semi_annual',
#     'high_activity_MSM_annual',
#     'high_activity_MSM_semi_annual',
#     "MSM_annual",
#     "MSM_semi_annual",
#     "high_activity_annual",
#     "high_activity_semi_annual",
#     "female_young_annual",
#     "female_young_semi_annual",
#     "remove_ltfu_10pct",
#     "remove_ltfu_20pct",
#     "remove_f_msm_10pct_and_msw_20pct_ltfu",
#     "1year_high_activity_semi_annual"
#   )
#
#   results <- list()
#   results[['BA']] <- list()
#   results[['SF']] <- list()
#
#   sf_year_idx <- years - 1942 + 1  # +1 to adjust for year 0, i.e. the initpop
#   ba_year_idx <- sf_year_idx + 1 # +1 for ba because we have an extra year of data for ba
#   # 1942 = gc_env$start.year - gc_env$cal.start
#   n_intervention_years <- max(years) - (gc_env$start.year + gc_env$cal.period)
#
#   # for (site in c("SF", "BA")) {
#     gc_env$model.end <- gc_env$model.end + n_intervention_years + 1
#     year_idx <- switch(site, SF = sf_year_idx, BA = ba_year_idx)
#
#     trace <- switch(
#       site,
#       BA = as.data.frame(readRDS(system.file("calibration_outcomes/BA_post_sample_strict_nhanes_targets.rds", package = "gcRegional"))$theta.list),
#       SF = as.data.frame(readRDS(system.file("calibration_outcomes/SF_post_sample_strict_nhanes_targets.rds", package = "gcRegional"))$theta.list)
#     )
#
#     for (i in 1:nrow(trace)) {
#       theta <- unlist(trace[i, ])
#       sols <- list()
#       results[[site]][[i]] <- list()
#       # simulate each intervention for a parameter vector
#       for (interv in interventions) {
#         results[[site]][[i]][[interv]] <- list()
#         e <- create_gcSim_param_env(theta, interv = interv)
#         run_gcSim_with_environment(e)
#         sols[[interv]] <- e$out.cpp$out
#         rm(e)
#       }
#       # now we compute the summary statistics for each simulation
#       for (interv in interventions) {
#         y.z.index <- c(1 + gc_env$y.index, 1 + gc_env$z.index)
#         s.y.z.index <- c(1 + gc_env$s.index, y.z.index)
#         prevalence <- rowSums(sols[[interv]][year_idx, y.z.index]) / rowSums(sols[[interv]][year_idx, s.y.z.index])
#         tests <- rowSums(sols[[interv]][year_idx, gc_env$scr.index]) - rowSums(sols[[interv]][min(year_idx)-1, gc_env$scr.index, drop=F])
#         results[[site]][[i]][[interv]][['prevalence']] <- prevalence
#         results[[site]][[i]][[interv]][['tests']] <- tests
#         if (interv != "base_case") {
#           # cases averted: how many cases occured in the base case but not in the intervention case?
#           base_case_cumulative_incidence <- rowSums(sols[["base_case"]][year_idx, gc_env$inc.index])
#           intervention_cumulative_incidence <- rowSums(sols[[interv]][year_idx, gc_env$inc.index])
#           cases_averted <- base_case_cumulative_incidence - intervention_cumulative_incidence
#           # number of tests necessary to avert one additional case
#           base_case_tests <- rowSums(sols[["base_case"]][year_idx, gc_env$scr.index]) - rowSums(sols[['base_case']][min(year_idx)-1, gc_env$scr.index, drop=F])
#           nnt <- (tests - base_case_tests) / (cases_averted)
#           results[[site]][[i]][[interv]][['cases_averted']] <- cases_averted
#           results[[site]][[i]][[interv]][['nnt']] <- nnt
#         }
#       }
#       if (i %% 100 == 0 && !missing(partial_results_path)) {
#         saveRDS(results, partial_results_path)
#       }
#     }
#   # }
#   return(results)
# }
