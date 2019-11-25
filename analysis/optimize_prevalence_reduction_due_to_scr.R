

library(gcRegional)

load_start('BA')
e <- create_gcSim_param_env(gc_env$theta)
gc_env$n_intervention_years <- max(nrow(e$screen) - (gc_env$cal.end), 0)
gc_env$intervention_rows <- (nrow(e$screen) - gc_env$n_intervention_years):nrow(e$screen)

#' Use a Vector of Percentages of Total Testing to Determine Prevalence
#' Reduction After 5 Years
#'
#' For each compartment of the population (1-28) we assign a screening rate to
#' be used in simulation for the 5 years after the data calibration period.
prevalence_after_screening <- function(scr = rep(1, 28), theta = gc_env$theta) {
  stopifnot(exists('site', envir = gc_env))
  e <- create_gcSim_param_env(theta)
  n_intervention_years <- max(nrow(e$screen) - (gc_env$cal.end), 0)
  intervention_rows <- (nrow(e$screen) - n_intervention_years):nrow(e$screen)
  for (row in intervention_rows) {
    e$screen[row, 1:28] <- scr
  }
  e$params$screen <- e$screen
  run_gcSim_with_environment(e)
  sol <- e$out.cpp$out
  last_year_prev <- sum(sol[nrow(sol), 1 + c(gc_env$y.index, gc_env$z.index)]) / sum(sol[nrow(sol), 1+c(gc_env$s.index, gc_env$y.index, gc_env$z.index)])
  return(last_year_prev)
}

# Get the Screenable Population Size in the Last Data Year
get_popsizes <- function(sol) {
  (sol[gc_env$cal.end, 1+gc_env$s.index] + sol[gc_env$cal.end, 1+gc_env$z.index])[1:28]
}


# Start with a log-fraction of screening that
# each population will receive, and back-calculate their
# screening rate.
convert_tau <- function(tau, popsizes, ntests) {
  # convert from log-scale to real > 0, and then
  # normalize so that each population gets a
  # percentage of the total amount of screening
  tau <- exp(tau)/sum(exp(tau))
  scr <- tau * ntests / popsizes
  return(scr)
}

transform_scr <- function(scr, popsizes) {
  # convert screening rates to fractions of screening
  # received
  tau <- scr * popsizes
  tau <- tau / sum(tau)
  tau <- log(tau)
  return(tau)
}


#' Optimize Screening
prevalence_reduction_via_screening <- function(tau, theta, increase = 0) {

  # construct simulation environment
  e <- create_gcSim_param_env(theta)

  # run simulation with theta in base_case
  run_gcSim_with_environment(e)
  sol <- e$out.cpp$out

  # calculate screenable popsize as symptomatic + asymptomatic
  popsizes <- get_popsizes(sol)

  # get screening rates as determined by theta
  theta_scr <- e$screen[gc_env$cal.end, 1:28]

  # get number of tests done in cal.end year
  ntests <- sum(theta_scr*popsizes)

  # get prevalence in last year of basecase simulation
  basecase_prev <- sum(sol[nrow(sol), 1 + c(gc_env$y.index, gc_env$z.index)]) / sum(sol[nrow(sol), 1+c(gc_env$s.index, gc_env$y.index, gc_env$z.index)])

  # convert tau from log-scale to >0 and normalize so that after the
  # cal.end year the number of tests is increased by approximately (increase)%.
  scr <- convert_tau(tau, popsizes, ntests * (1+increase))

  # update screening matrix
  for (row in gc_env$intervention_rows) {
    e$screen[row, 1:28] <- scr
  }
  e$params$screen <- e$screen

  # re-run model
  run_gcSim_with_environment(e)
  sol <- e$out.cpp$out

  # get prevalence in last year after changing prevalence
  scr_prev <- sum(sol[nrow(sol), 1 + c(gc_env$y.index, gc_env$z.index)]) / sum(sol[nrow(sol), 1+c(gc_env$s.index, gc_env$y.index, gc_env$z.index)])

  # return:
  #  - relative reduction in prevalence
  return(list(
    prevalence_relative_reduction = (basecase_prev - scr_prev) / basecase_prev,
    scr = scr,
    relative_scr = scr / theta_scr
  ))
}

optim_fun_0 <- function(tau, increase = 0) {
  # We want exp(tau) to reflect the proportion of tests
  # given to each demographic, and the total = 1.
  # Even if exp(tau) != 1, convert_tau will normalize it,
  # and so the correct proportions will still be used.
  # The point of the penalty is to make sure that if there
  # is a single optimum when sum(exp(tau))==1 that this optimum
  # doesn't turn into a ray in the optimization space.
  penalty <- (1-sum(exp(tau)))^2
  out <- prevalence_reduction_via_screening(tau, theta = theta, increase = 0)
  val <- out$prevalence_relative_reduction - penalty
  print(val)
  return(val)
}


theta <- gc_env$theta

# construct simulation environment
e <- create_gcSim_param_env(theta)

# run simulation with theta in base_case
run_gcSim_with_environment(e)
sol <- e$out.cpp$out

# calculate screenable popsize as symptomatic + asymptomatic
popsizes <- get_popsizes(sol)

# get screening rates as determined by theta
theta_scr <- e$screen[gc_env$cal.end, 1:28]

# get number of tests done in cal.end year
ntests <- sum(theta_scr*popsizes)

tau <- theta_scr*popsizes
tau <- tau / sum(tau)
tau <- log(tau)

out <- optim(par = tau, fn = optim_fun_0, control = c(fnscale = -1))
