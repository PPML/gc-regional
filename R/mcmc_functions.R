##############################################
###   additional functions used for MCMC  ####
###   and other miscellaneous functions   ####
##############################################

#calculate likelihood of given iteration and update parameter set used
likelihood <- function(theta){
  if(check_variation('imis')) {
    return(apply(theta, 1, function(x) exp(model_epi_loglik(x)[['ll']])))
  }
  if (check_variation('national_posterior_medians')) {
    theta <- replace_national_parameter_medians(theta)
  }
  if (check_variation('no_differential_reporting')) {
    theta[c("logit.risk.rep.symp.m",  "logit.risk.rep.symp.m1", "logit.risk.rep.symp.f",  "logit.risk.rep.symp.f1")] <- logit(1)
  }
  Res <- suppressWarnings(model_epi_loglik(theta)) #run the model
  LogLL <- Res["ll"] #Model is returning the Log likelihood
  theta.last<-theta #update theta
  gc_assign(theta.last)
  out <- list( LogLL=LogLL)
  return(out)
}

#calculate posterior likelihood
#' @export
dLogPosterior <- function(theta) {
  tryCatch({
    if (is.null(names(theta))) {
      names(theta) <- names(gc_env$theta)
    }

    log.prior <- prior(theta)
    log.likelihood <- unlist(likelihood(theta)$LogLL)
    log.posterior <- log.prior + log.likelihood

    if (is.na(log.posterior)||is.nan(log.posterior)) log.posterior <- -1e32
    return(as.numeric(log.posterior))
  },
  error = function(x) -1e32)
}



#run adaptive Metroposlis Hastings MCMC algorithm (adapted from fitR package)
#see mcmcMH in fitR package for details
#' @export
my_mcmc<-function(target, init.theta, proposal.sd, n.iterations, print.info.every = n.iterations / 100) {
  out_all <- NULL
  gc_assign(out_all)
  if (check_variation('no_adaptive')) {
    adapt.size.start = n.iterations
    adapt.shape.start = n.iterations
  } else {
    adapt.size.start = 1000
    adapt.shape.start = 3000
  }
  trace <- myMCMC(target = target,
                  init.theta = init.theta,
                  proposal.sd = proposal.sd,
                  n.iterations = n.iterations,
                  adapt.size.start = adapt.size.start,
                  adapt.shape.start = adapt.shape.start,
                  adapt.size.cooling=0.999,
                  print.info.every = print.info.every,
                  verbose=FALSE)
  return(trace)
}

# function to estimate parameters for beta distribution, where pred.prop =
# probability, pop = population size, and scale is a scaling factor that
# determines variance
estBetaParams <- function(pred.prob, pop, scale) {
  alpha <- pred.prob*pop/scale
  beta <- (1-pred.prob)*pop/scale
  return(beta.params = list(alpha = alpha, beta = beta))
}

estBetaParamsVar <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  list(alpha = alpha, beta = beta)
}

estBetaMeanVar <- function(shape1, shape2) {
  mean = shape1 / (shape1 + shape2)
  var = shape1 * shape2 / ( (shape1 + shape2)^2 * (shape1 + shape2 + 1))
  return(c(mean = mean, var = var))
}


logit<-function(x) {log(x/(1-x))}

ilogit <-function(x) {1/(1+exp(-x))}


# function to calculate annual_rates from cumulative output, when model is
# outputting ANNUAL values -- would need to modify if outputting weekly or
# monthly values
annual_rates <- function(x){
  annual <-rbind(rep(0,ncol(x)), apply(x,2,diff))
}

