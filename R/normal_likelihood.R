big_sd_normal_likelihoods <- function(pred, N = 5, log = T) {
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

  likelihood <- list()

  max_p_symp_sd <- .6/N
  likelihood[['p.symp']] <-
    if(site == "SF") {
      dnorm(
        x = as.numeric(unlist(pred[["prev"]][["fit.symp"]]))[2:length(p.symp.ssun)],
        mean = p.symp.ssun[2:length(p.symp.ssun)],
        sd = max_p_symp_sd*N,
        log = log
      )
    } else {
      dnorm(
        x = as.numeric(unlist(pred[["prev"]][["fit.symp"]])),
        mean = p.symp.ssun,
        sd = max_p_symp_sd*N,
        log = log
      )
    }

  max_diag_rate_sd <- max(diag.rate.sd)*1.2
  likelihood[['diag.m']] <-
    dnorm(
      x = as.numeric(unlist(pred[["prev"]][["fit.diag.rate"]]))[1:(length(diag.rate)/2)],
      mean = (diag.rate)[1:(length(diag.rate)/2)],
      sd = max_diag_rate_sd*N,
      log = log
    )

  likelihood[['diag.f']] <-
    dnorm(
      x = as.numeric(unlist(pred[["prev"]][["fit.diag.rate"]]))[((length(diag.rate)/2)+1):length(diag.rate)],
      mean = (diag.rate)[((length(diag.rate)/2)+1):length(diag.rate)],
      sd = max_diag_rate_sd*N,
      log = log
    )

  max_rr_diag_subpop_sd <- max(rr.diag.subpop.sd)*1.2
  likelihood[['diag.rr']] <-
    dnorm(
      x = as.numeric(unlist(pred[["prev"]][["diag.rr"]])),
      mean = rr.diag.subpop,
      sd = max_rr_diag_subpop_sd*N,
      log = log
    )

  # max_p_msm <- max(sqrt(var.p.msm.ssun))
  max_p_msm <- .25/N
  likelihood[['p.msm']] <-
    dnorm(
      x = as.numeric(unlist(pred[["prev"]][["p.diag.msm"]])),
      mean = p.msm.ssun,
      sd = max_p_msm*N,
      log = log
    )

  max_age_sd <- max(sqrt(age.dist.dat$var))
  likelihood[['age']] <-
    dnorm(
      x = pred[["age.dist.all"]][1:4],
      mean = (age.dist.dat$p.same.age)[1:4],
      sd = max_age_sd*N,
      log = log
    )
  return(likelihood)
}
