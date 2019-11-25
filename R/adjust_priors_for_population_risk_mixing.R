adjust_priors_for_population_risk_mixing <- function() {
  stopifnot(exists('priors', envir = gc_env))
  gc_env$priors <- rbind.data.frame(gc_env$priors,
    data.frame(
      parameter = 'epsilon',
      distribution = 'beta',
      param1 = 1.1,
      param2 = 1.1,
      mean = .5,
      sd = .39,
      sd.prior = .59,
      transformation = 'logit',
      mean.transf = 0,
      ucl.transf = 3.411,
      sd.transf.1 = -3.411,
      sd.transf = -3.411,
      sd1 = -0.45
    ))

  gc_env$priors[gc_env$priors$parameter %in% c('epsilon.1', 'epsilon.2', 'epsilon.3', 'epsilon.4'), 'param1'] <- 100
  gc_env$priors[gc_env$priors$parameter %in% c('epsilon.1', 'epsilon.2', 'epsilon.3', 'epsilon.4'), 'param2'] <- 100

}
