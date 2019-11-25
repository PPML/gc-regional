# Optimize Prevalence Reduction Via Screening

# Given theta, a scaling factor to apply to the current total level of screening,
# how can we determine the most effective redistribution of screening?

# We will have two ways of doing this:
# Either we will have no identification of the High Activity groups, or we will
# apply the adjusted screening rates for the High-Activity group to 50% of the
# High Activity group and 10% of the Low Activity group in the same demographic
# group.




# Helper Functions --------------------------------------------------------

# Get the Screenable Population Size in the Last Data Year
get_popsizes <- function(sol) {
  (sol[gc_env$cal.end, 1+gc_env$s.index] + sol[gc_env$cal.end, 1+gc_env$z.index])[1:28]
}


transform_scr <- function(scr, popsizes) {
  # convert screening rates to fractions of screening
  # received
  tau <- scr * popsizes
  tau <- tau / sum(tau)
  tau <- log(tau)
  return(tau)
}



# Optimize Screening ------------------------------------------------------

optimize_screening <-
  function(theta,
           identifiability = c('none',
                               'partial',
                               'complete'),
           increase = 0) {

    # pre-flight inspection
    #   - a site needs to have been load_start()'ed
    if (! exists('site', envir = gc_env)) stop("Site-specific data must be loaded.")

    #   - identifiability must be specified properly
    if (! identifiability %in% c('none',
                                 'partial',
                                 'complete') &&
        length(identifiability) == 1) {
      stop("identafiability must be one of none, partial, or complete.")
    }

    # construct the simulation environment containing the screening matrix from theta
    e <- create_gcSim_param_env(theta)

    # get the time-frame for our intervention
    n_intervention_years <- max(nrow(e$screen) - (gc_env$cal.end), 0)
    intervention_rows <- (nrow(e$screen) - n_intervention_years):nrow(e$screen)

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
    basecase_prev <- sum(sol[nrow(sol), 1 + c(gc_env$y.index, gc_env$z.index)]) /
      sum(sol[nrow(sol), 1+c(gc_env$s.index, gc_env$y.index, gc_env$z.index)])

    # Construct tau from theta_scr
    tau <- theta_scr*popsizes
    tau <- tau / sum(tau)
    tau <- log(tau)

    # Start with a log-fraction of screening that
    # each population will receive, and back-calculate their
    # screening rate.
    convert_tau <- function(tau, popsizes, ntests) {

      # convert from log-scale to real > 0, and then
      # normalize so that each population gets a
      # percentage of the total amount of screening
      tau <- exp(tau)/sum(exp(tau))
      scr <- tau * ntests / popsizes

      # apply identifiability
      switch(identifiability,
             none = {
               # replace high activity with low activity rates
               scr[0:13*2+2] <- scr[0:13*2+1]
             },
             partial = {
               # low activity gets 10% of screening intended for high activity
               scr[0:13*2+1] <- (.9 * scr[0:13*2+1]) + (.1 * scr[0:13*2+2])
               # high activity gets half of the intended low activity rate and
               # half of the intended high activity rate
               scr[0:13*2+2] <- (.5 * scr[0:13*2+1]) + (.5 * scr[0:13*2+2])
             },
             complete = {} # do nothing
      )
      return(scr)
    }

    # simulate and output the prevalence reduction from an altered and transformed screening vector
    prevalence_after_screening <- function(tau) {

      # convert tau from log-scale to >0 and normalize so that after the
      # cal.end year the number of tests is increased by approximately (increase)%.
      scr <- convert_tau(tau, popsizes, ntests * (1+increase))

      # update screening matrix
      for (row in intervention_rows) {
        e$screen[row, 1:28] <- scr
      }
      e$params$screen <- e$screen

      # run model
      run_gcSim_with_environment(e)
      sol <- e$out.cpp$out

      # get prevalence in last year after changing prevalence
      scr_prev <- sum(sol[nrow(sol), 1 + c(gc_env$y.index, gc_env$z.index)]) /
        sum(sol[nrow(sol), 1+c(gc_env$s.index, gc_env$y.index, gc_env$z.index)])

      return(list(
        prevalence_relative_reduction = (basecase_prev - scr_prev) / basecase_prev,
        scr = scr,
        relative_scr = scr / theta_scr
      ))
    }

    # Reward Function for Optimization
    optim_fun_0 <- function(tau) {
      # We want exp(tau) to reflect the proportion of tests
      # given to each demographic, and the total = 1.
      # Even if exp(tau) != 1, convert_tau will normalize it,
      # and so the correct proportions will still be used.
      # The point of the penalty is to make sure that if there
      # is a single optimum when sum(exp(tau))==1 that this optimum
      # doesn't turn into a ray in the optimization space.
      penalty <- (1 + increase -sum(exp(tau)))^2
      out <- prevalence_after_screening(tau)
      val <- out$prevalence_relative_reduction - penalty
      print(val)
      return(val)
    }

    # run optimization
    optim(par = tau, fn = optim_fun_0, control = c(fnscale = -1))
  }
