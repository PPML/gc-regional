##################################################################
### function to use model outputs to calculate log likelihoods ###
##################################################################
model_epi_loglik <- function(theta, debug=F) {
  # retrieve environment variables
  site <- gc_env$site
  p.symp.ssun <- gc_env$p.symp.ssun
  diag.rate <- gc_env$diag.rate
  rr.diag.subpop <- gc_env$rr.diag.subpop
  rr.diag.subpop.sd <- gc_env$rr.diag.subpop.sd
  p.msm.ssun <- gc_env$p.msm.ssun
  age.dist.dat <- gc_env$age.dist.dat
  var.symp.ssun <- gc_env$var.symp.ssun
  var.p.msm.ssun <- gc_env$var.p.msm.ssun
  diag.rate.sd <- gc_env$diag.rate.sd
  var.prev.nhanes <- gc_env$var.prev.nhanes
  prev.nhanes <- gc_env$prev.nhanes
  nhanes.updated.dat <- gc_env$nhanes.updated.dat

  if (widen_likelihood()) {
    N <- 2
    var.prev.nhanes <- var.prev.nhanes*N
    diag.rate <- diag.rate*sqrt(N)
    age.dist.dat$var <- age.dist.dat$var*N
    var.symp.ssun <- var.symp.ssun*N
    var.p.msm.ssun <- var.p.msm.ssun*N
  }

  if (check_variation('wide_normal_likelihoods')) {
    pred <- prediction_epi(theta)
    lik <- big_sd_normal_likelihoods(pred)
    ll <- sum(sapply(lik, sum))
    return(list(ll = ll))
  }

  if (check_variation('wide_likelihood')) {
    pred <- prediction_epi(theta)
    lik <- wide_likelihood(pred)
    ll <- sum(sapply(lik, sum))
    return(list(ll = ll))
  }

  if (check_variation('population_denominators')) {
    pred <- prediction_epi(theta)
    lik <- population_weighted_likelihood(pred)
    ll <- sum(sapply(lik, sum))
    return(c(unlist(pred), ll=ll, unlist(lik)))
  }

  pred <- prediction_epi(theta) # Run the model

  ### calculate likelihoods
  #for beta parameter estimates, using model estimate as mean, data estimates for variance
  #beta distn parameter estimates for age assortativity
  beta.params.age <- estBetaParamsVar(pred[["age.dist.all"]],age.dist.dat$var)

  # We are no longer using the original nhanes estimates data. Updated nhanes data is
  # being used as a lower bound on the prevalence rates in each compartment.
  # beta.params.prev <- estBetaParamsVar(as.numeric(unlist(pred[["prev"]][["fit.prev.extra"]])),
  #                  var.prev.nhanes) #fit to extra prev cats, include MSM in estimate

  #beta distn parameter estimates for proportion symptomatic reported cases
  beta.params.symp <- estBetaParamsVar(as.numeric(unlist(pred[["prev"]][["fit.symp"]])), var.symp.ssun)

  #beta distn parameter estimates for reported cases
  beta.params.diag <- estBetaParamsVar(as.numeric(unlist(pred[["prev"]][["fit.diag.rate"]])), diag.rate.sd^2 )

  #beta distn parameter estimates for proportion of males cases who are MSM
  beta.params.pmsm <- estBetaParamsVar(as.numeric(unlist(pred[["prev"]][["p.diag.msm"]])), var.p.msm.ssun)

  #log likelihood fitting to overall pop estimates of age assortativity only
  ll.age <- sum(dbeta(x=age.dist.dat$p.same.age[1:4],beta.params.age$alpha[1:4], beta.params.age$beta[1:4], log=TRUE))/2

  # log likelihood symptomatic
  ll.symp <- ifelse(
    site != "SF",
    sum(dbeta(
      x = p.symp.ssun,
      beta.params.symp$alpha,
      beta.params.symp$beta,
      log = TRUE
    )),
    sum(dbeta(
      x = p.symp.ssun[2:length(p.symp.ssun)],
      beta.params.symp$alpha[2:length(p.symp.ssun)],
      beta.params.symp$beta[2:length(p.symp.ssun)],
      log = TRUE
    )) / 3)

  #log likelihood reported case rate male
  ll.diag.m <- sum(dbeta(x=diag.rate[1:(length(diag.rate)/2)],
  beta.params.diag$alpha[1:(length(beta.params.diag$alpha)/2)],beta.params.diag$beta[1:(length(beta.params.diag$beta)/2)],
  log=TRUE))

  #log likelihood reported case rate female
  ll.diag.f <-
  sum(dbeta(x=diag.rate[((length(diag.rate)/2)+1):length(diag.rate)],
  beta.params.diag$alpha[((length(beta.params.diag$alpha)/2)+1):length(beta.params.diag$alpha)],
  beta.params.diag$beta[((length(beta.params.diag$beta)/2)+1):length(beta.params.diag$beta)],
  log=TRUE))

  #log likelihood reported case rr by subpop
  ll.diag.subpop  <- sum(dnorm(x=rr.diag.subpop, mean=as.numeric(unlist(pred[["prev"]][["diag.rr"]])),sd=rr.diag.subpop.sd*2, log=TRUE))/5

  # Log Likelihood Prevalence
  if (exists('gc_testing', envir = .GlobalEnv) && 'nhanes_penalty' %in% gc_testing) {
  # Nhanes penalty term:
  # This essentially makes being less than the NHANES prevalence estimates
  # a penalty orders of magnitude higher than any other factor in the log posterior.
    nhanes_diff <- pred[['prev']][['msm_redistributed_prevalence_rates']] - nhanes.updated.dat$prev
    ll.nhanes.prev <- -10^3*sum(abs(nhanes_diff[which(nhanes_diff < 0)]))
  } else if (exists('gc_testing', envir = .GlobalEnv) && 'nhanes_lower_uniform' %in% gc_testing) {
    ll.nhanes.prev <- ifelse(
      all(
        (nhanes.updated.dat$prev - nhanes.updated.dat$prev_std) <=
          pred[['prev']][['msm_redistributed_prevalence_rates']]),
      0, -Inf)
  } else if (exists('gc_testing', envir = .GlobalEnv) && 'no_nhanes_limit' %in% gc_testing) {
    ll.nhanes.prev <- 0
  } else {
  # We require that the prevalence rates in categories that can be compared to the NHANES data
  # be at least as high as the NHANES estimate. log(0)=-Inf and log(1)=0, so when any
  # category has prevalence < NHANES, the ll.nhanes.prev = -Inf, and when all categories have
  # prevalence >= NHANES, ll.nhanes.prev = 0.
    ll.nhanes.prev <- ifelse(
      all(nhanes.updated.dat$prev <= pred[['prev']][['msm_redistributed_prevalence_rates']]),
      0, -Inf)
  }



  #log likelihood pMSM
  ll.msm <-sum(dbeta(x=p.msm.ssun,beta.params.pmsm$alpha, beta.params.pmsm$beta, log=TRUE))/2

  # log likelihood of msm prevalence
  # for now these are hardcoded here, but this needs to changed so that they're
  # explicitly declared in a parameter table.
  ll.msm.prev <- 0
  if (gc_env$site == "SF") {
    if (exists('gc_testing', envir = .GlobalEnv)) {
      if ('widen_msm_prevalence' %in% gc_testing) {
        msm_prev_alpha <- 5
        msm_prev_beta <- 78.33
        ll.msm.prev <- dbeta(x=pred[['prev']][['fit.prev.extra']][['prev.msm']], msm_prev_alpha, msm_prev_beta, log=T)
      } else if('no_msm_prevalence' %in% gc_testing) {
        ll.msm.prev <- 0
      } else {
        msm_prev_alpha <- 89.744
        msm_prev_beta <- 1405.989
        ll.msm.prev <- dbeta(x=pred[['prev']][['fit.prev.extra']][['prev.msm']], msm_prev_alpha, msm_prev_beta, log=T)
      }
    }

  } else ll.msm.prev <- 0

  #sum the log likihoods
  ll.terms <-
    c(
      ll.age = ll.age,
      ll.diag.f = ll.diag.f,
      ll.diag.m = ll.diag.m,
      ll.diag.subpop = ll.diag.subpop,
      ll.msm = ll.msm,
      ll.symp = ll.symp,
      ll.nhanes.prev = ll.nhanes.prev,
      ll.msm.prev = ll.msm.prev
    )
  if (debug) {
    if (any(is.na(ll.terms))) {
      warning("log likelihood contained NAs")
      print(ll.terms)
    }
  }
  ll <- sum(ll.terms)
  # ll<-sum(ll.age, ll.diag.f, ll.diag.m, ll.diag.subpop, ll.msm, ll.symp, ll.nhanes.prev)

  #to prevent errors
  ll[is.na(ll)]<-(-1e20)

  #return model outputs
  c(unlist(pred), ll.age=ll.age, ll.symp=ll.symp, ll.diag.m=ll.diag.m,
  ll.diag.f=ll.diag.f, ll.diag.subpop=ll.diag.subpop, ll.msm=ll.msm, ll.nhanes.prev=ll.nhanes.prev, ll=ll)
}
