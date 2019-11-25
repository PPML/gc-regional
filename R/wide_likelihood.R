wide_likelihood <- function(pred, N = 5, log = T) {
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

  likelihood <- list()

  likelihood[['p.symp']] <-
    if(site == "SF") {
      dbeta(
        x = p.symp.ssun[2:length(p.symp.ssun)],
        shape1 = (as.numeric(unlist(pred[["prev"]][["fit.symp"]]))[2:length(p.symp.ssun)])*100,
        shape2 = (1-as.numeric(unlist(pred[["prev"]][["fit.symp"]]))[2:length(p.symp.ssun)])*100,
        log = log
      )
    } else {
      dbeta(
        x = p.symp.ssun,
        shape1 = as.numeric(unlist(pred[["prev"]][["fit.symp"]]))*100,
        shape2 = (1-as.numeric(unlist(pred[["prev"]][["fit.symp"]])))*100,
        log = log
      )
    }


  likelihood[['p.msm']] <-
    dbeta(
      x = p.msm.ssun,
      shape1 = as.numeric(unlist(pred[["prev"]][["p.diag.msm"]]))*100,
      shape2 = (1-as.numeric(unlist(pred[["prev"]][["p.diag.msm"]])))*100,
      log = log)

  likelihood[['age']] <-
    dbeta(
      x = (age.dist.dat$p.same.age)[1:4],
      shape1 = pred[["age.dist.all"]][1:4]*100,
      shape2 = (1-pred[["age.dist.all"]][1:4])*100,
      log = log)

  nhanes_dbeta <- function(N) {
    dbeta(
      x = nhanes.updated.dat$prev,
      shape1 = pred$prev[['msm_redistributed_prevalence_rates']]*N,
      shape2 = (1-pred$prev[['msm_redistributed_prevalence_rates']])*N,
      log = log
    ) 
  }

  # In the NHANES targets, we exclude the 15th likelihood component
  # because it's a 0% prevalence estimate for the 25-39 female 
  # hispanic population
  if (check_variation('use_nhanes_50')) {
    likelihood[['nhanes_prev']] <- nhanes_dbeta(50)[-15]
  } 
  if (check_variation('use_nhanes_100')) {
    likelihood[['nhanes_prev']] <- nhanes_dbeta(100)[-15]
  } 
  if (check_variation('just_male_nhanes_targets')) {
    male_indices <- nhanes.updated.dat$sex == 'M'
    likelihood[['nhanes_prev']] <- nhanes_dbeta(50)[male_indices]
  }
  if (check_variation('just_female_nhanes_targets')) {
    female_indices <- which(nhanes.updated.dat$sex == 'W')
    female_indices <- setdiff(female_indices, 15)
    likelihood[['nhanes_prev']] <- nhanes_dbeta(50)[female_indices]
  }

  if (check_variation('use_case_rates')) {
    likelihood[['case_rates']] <-
      dbeta(
        x = unlist(c(diag.subpop.rate)),
        shape1 = pred$prev[['fit.diag.subpop']]*100,
        shape2 = (1-pred$prev[['fit.diag.subpop']])*100,
        log = log
      )
  } else {
    likelihood[['diag.m']] <-
      dbeta(
        x = (diag.rate)[1:(length(diag.rate)/2)],
        shape1 = as.numeric(unlist(pred[["prev"]][["fit.diag.rate"]]))[1:(length(diag.rate)/2)]*100,
        shape2 = (1-as.numeric(unlist(pred[["prev"]][["fit.diag.rate"]]))[1:(length(diag.rate)/2)])*100,
        log = log)

    likelihood[['diag.f']] <-
      dbeta(
        x = (diag.rate)[((length(diag.rate)/2)+1):length(diag.rate)],
        shape1 = as.numeric(unlist(pred[["prev"]][["fit.diag.rate"]]))[((length(diag.rate)/2)+1):length(diag.rate)]*100,
        shape2 = (1-as.numeric(unlist(pred[["prev"]][["fit.diag.rate"]]))[((length(diag.rate)/2)+1):length(diag.rate)])*100,
        log = log)

    likelihood[['diag.rr']] <-
      dnorm(
        x = rr.diag.subpop,
        mean = as.numeric(unlist(pred[["prev"]][["diag.rr"]])),
        sd = rr.diag.subpop.sd * 2,
        log = TRUE
      ) / 5

  }

  return(likelihood)
}
