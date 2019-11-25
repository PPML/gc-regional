#' Population Weighted Likelihood Function
population_weighted_likelihood <- function(pred) {
  # get target data from gc_env ---------------------------------------------
  p.symp.ssun <- gc_env$p.symp.ssun
  diag.rate <- gc_env$diag.rate
  rr.diag.subpop <- gc_env$rr.diag.subpop
  rr.diag.subpop.sd <- gc_env$rr.diag.subpop.sd
  p.msm.ssun <- gc_env$p.msm.ssun
  age.dist.dat <- gc_env$age.dist.dat
  var.symp.ssun <- gc_env$var.symp.ssun
  var.p.msm.ssun <- gc_env$var.p.msm.ssun
  diag.rate.sd <- gc_env$diag.rate.sd
  site <- gc_env$site
  nhanes.updated.dat <- gc_env$nhanes.updated.dat
  diag.subpop.rate <- gc_env$diag.subpop.rate
  p.symp.popsizes <- gc_env$p.symp.popsizes
  observed_proportion_msm_cases <- gc_env$observed_proportion_msm_cases
  n.i <- gc_env$n.i
  m4 <- gc_env$m4

  if (exists('scalars', envir = gc_env)) scalars <- gc_env$scalars
  else scalars <- c(p.symp = 1, p.msm = 1, age = 1, nhanes_prev = 1, case_rates = 1)

  # get scalars based on named lookup in the scalars parameter vector
  get_scalar <- function(str) {
    if (str %in% names(scalars)) return(scalars[[str]])
    else return(1)
  }

  # set up likelihood
  likelihood <- list()

  # p.symp - proportion symptomatic --------------------------------------------------
  predicted_p_symp <- as.numeric(unlist(pred[["prev"]][["fit.symp"]]))

  likelihood[['p.symp']] <-
    dbeta(
      x = p.symp.ssun,
      predicted_p_symp * p.symp.popsizes * get_scalar("p.symp"),
      (1-predicted_p_symp) * p.symp.popsizes * get_scalar("p.symp"),
      log = TRUE
    )
  if (site == 'SF') likelihood[['p.symp']] <- likelihood[['p.symp']][2:length(likelihood[['p.symp']])]

  # p.msm - proportion of cases among msm -------------------------------------------
  predicted_p_diag_msm <- as.numeric(unlist(pred[["prev"]][["p.diag.msm"]]))

  likelihood[['p.msm']] <- dbeta(
    x = p.msm.ssun,
    predicted_p_diag_msm * observed_proportion_msm_cases,
    (1-predicted_p_diag_msm) * observed_proportion_msm_cases,
    log = T
  )

  # age - age assortativity -------------------------------------------------------
  predicted_age_dist <- pred[["age.dist.all"]]
  likelihood[['age']] <-
    dbeta(
      x = age.dist.dat$p.same.age,
      shape1 = predicted_age_dist*age.dist.dat$N,
      shape2 = (1-predicted_age_dist)*age.dist.dat$N,
      log=T
    )
  # beta.params.age <- estBetaParamsVar(pred[["age.dist.all"]],age.dist.dat$var)
  # likelihood[['age']] <- dbeta(x=age.dist.dat$p.same.age[1:4],beta.params.age$alpha[1:4 * get_scalar('age')], beta.params.age$beta[1:4] * get_scalar('age'), log=TRUE)

  # nhanes_prev - nhanes prevalence estimates ---------------------------------------------
  if (! check_variation('dont_use_nhanes')) {
    likelihood[['nhanes_prev']] <- dbeta(
      x = nhanes.updated.dat$prev,
      shape1 = pred$prev[['msm_redistributed_prevalence_rates']]*nhanes.updated.dat$prev_denom * get_scalar('nhanes_prev'),
      shape2 = (1-pred$prev[['msm_redistributed_prevalence_rates']])*nhanes.updated.dat$prev_denom * get_scalar('nhanes_prev'),
      log = T
    )[-15] # exclude the 0 prevalence estimate for 25-39 hispanic W

    if (check_variation('only_use_nhanes_women')) {
      likelihood[['nhanes_prev']] <- likelihood[['nhanes_prev']][setdiff(which(gc_env$nhanes.updated.dat$sex == 'W'), 15)]
    }
  }


  # case_rates --------------------------------------------------------------
  likelihood[['case_rates']] <-
    dbeta(
      x = unlist(c(diag.subpop.rate)),
      shape1 = pred$prev[['fit.diag.subpop']]*unlist(gc_env$diag.subpop.rate.denom) * get_scalar('case_rates'),
      shape2 = (1-pred$prev[['fit.diag.subpop']])*unlist(gc_env$diag.subpop.rate.denom) * get_scalar('case_rates'),
      log = T
    )

  if (check_variation('msm_incidence_likelihood')) {
    nyrs <- length(pred$prev$inc)/14
    likelihood[['msm_incidence']] <- dgamma(
      x = mean(pred$prev$inc[1:(nyrs*2) + 6*nyrs]*100/sum(n.i[m4])),
      shape = 3.358014,
      rate = 0.389351,
      log = T
    ) * 10
  }

  return(likelihood)
}
