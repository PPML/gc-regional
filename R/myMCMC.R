#'  Metropolis Hastings MCMC
#'
#'  This is a slight adaptation of mcmcMH from fitR package
#'  prevents 'covariance matrix not positive definite' error
#'  see fitR package for documentation.
#'
#'  Original Function located here:
#'  https://sbfnk.github.io/fitR/reference/mcmcMH.html
#'
#' References for Adaptive MCMC:
#'
#' https://people.eecs.berkeley.edu/~jordan/sail/readings/andrieu-thoms.pdf
#' http://probability.ca/jeff/ftpdir/adaptex.pdf
#' https://projecteuclid.org/download/pdf_1/euclid.bj/1080222083

myMCMC<-function (target, init.theta, proposal.sd = NULL, n.iterations,
                  covmat = NULL, limits = list(lower = NULL, upper = NULL),
                  adapt.size.start = NULL, adapt.size.cooling = 0.99, adapt.shape.start = NULL,
                  print.info.every = n.iterations/100, verbose = FALSE, max.scaling.sd = 50,
                  acceptance.rate.weight = NULL, acceptance.window = NULL)
{
  cal.period<- gc_env$cal.period
  cal.start <- gc_env$cal.start
  m1 <- gc_env$m1
  m2 <- gc_env$m2
  m3<- gc_env$m3

  theta.current <- init.theta
  theta.propose <- init.theta
  covmat.proposal <- covmat
  lower.proposal <- limits$lower
  upper.proposal <- limits$upper
  theta.names <- names(init.theta)
  if (!is.null(proposal.sd) && is.null(names(proposal.sd))) {
    names(proposal.sd) <- theta.names
  }
  if (is.null(covmat.proposal)) {
    if (is.null(proposal.sd)) {
      proposal.sd <- init.theta/10
    }
    covmat.proposal <- matrix(diag(proposal.sd[theta.names]^2,
                                   nrow = length(theta.names)), nrow = length(theta.names),
                              dimnames = list(theta.names, theta.names))
  }
  else {
    covmat.proposal <- covmat.proposal[theta.names, theta.names]
  }
  if (is.null(lower.proposal)) {
    lower.proposal <- init.theta
    lower.proposal[] <- -Inf
  }
  else {
    lower.proposal <- lower.proposal[theta.names]
  }
  if (is.null(upper.proposal)) {
    upper.proposal <- init.theta
    upper.proposal[] <- Inf
  }
  else {
    upper.proposal <- upper.proposal[theta.names]
  }
  covmat.proposal.init <- covmat.proposal
  adapting.size <- FALSE
  adapting.shape <- FALSE
  theta.estimated.names <- names(which(diag(covmat.proposal) >
                                         0))
  target.theta.current <- target(theta.current)

  if (class(target.theta.current) == "numeric") {
    suppressWarnings(target.theta.current$log.density <- target.theta.current)
    suppressWarnings(target.theta.current$trace <- theta.current)
  }
  if (!is.null(print.info.every)) {
    message("Init: ", fitR::printNamedVector(theta.current[theta.estimated.names]),
            ", target: ", target.theta.current[["log.density"]])
  }
  trace <- data.frame(t(target.theta.current[["trace"]]))
  acceptance.rate <- 0
  if (!is.null(acceptance.window)) {
    acceptance.series <- c()
  }
  scaling.sd <- 1
  scaling.multiplier <- 1
  covmat.empirical <- covmat.proposal
  covmat.empirical[, ] <- 0
  theta.mean <- theta.current
  if (is.null(print.info.every)) {
    print.info.every <- n.iterations + 1
  }
  start_iteration_time <- Sys.time()
  for (i.iteration in seq_len(n.iterations)) {
    if (!is.null(adapt.size.start) && i.iteration >= adapt.size.start &&
        (is.null(adapt.shape.start) || acceptance.rate *
         i.iteration < adapt.shape.start)) {
      if (!adapting.size) {
        message("\n---> Start adapting size of covariance matrix")
        adapting.size <- TRUE
      }
      scaling.multiplier <- exp(adapt.size.cooling^(i.iteration -
                                                      adapt.size.start) * (acceptance.rate - 0.234))
      scaling.sd <- scaling.sd * scaling.multiplier
      scaling.sd <- min(c(scaling.sd, max.scaling.sd))
      covmat.proposal.new <- scaling.sd^2 * covmat.proposal.init
      #print(covmat.proposal)
      if (!(any(diag(covmat.proposal.new)[theta.estimated.names] <
                .Machine$double.eps))) {
        covmat.proposal <- covmat.proposal.new
      }
    }
    else if (!is.null(adapt.shape.start) && acceptance.rate *
             i.iteration >= adapt.shape.start) {
      if (!adapting.shape) {
        message("\n---> Start adapting shape of covariance matrix")
        adapting.shape <- TRUE
      }
      scaling.sd <- 2.38/sqrt(length(theta.estimated.names))
      covmat.proposal <- scaling.sd^2 * covmat.empirical
    }
    if (i.iteration%%ceiling(print.info.every) == 0) {
      state.mcmc <- trace[nrow(trace), ]
      message("Iteration: ", i.iteration, "/", n.iterations,
              ", acceptance rate: ", sprintf("%.3f", acceptance.rate),
              appendLF = FALSE)
      if (!is.null(adapt.size.start) || !is.null(adapt.shape.start)) {
        message(", scaling.sd: ", sprintf("%.3f", scaling.sd),
                ", scaling.multiplier: ", sprintf("%.3f",
                                                  scaling.multiplier), appendLF = FALSE)
      }
      # message(", state: ", fitR::printNamedVector(state.mcmc))
    }
    if (any(diag(covmat.proposal)[theta.estimated.names] <
            .Machine$double.eps)) {
      print(covmat.proposal[theta.estimated.names, theta.estimated.names])
      stop("non-positive definite covmat", call. = FALSE)
    }
    if (length(theta.estimated.names) > 0) {
      if(det(covmat.proposal)==0) {
        covmat.proposal <- covmat.proposal[theta.estimated.names,theta.estimated.names]+diag(ncol(covmat.proposal))*0.01 ### change made here
      }
      theta.propose[theta.estimated.names] <- as.vector(tmvtnorm::rtmvnorm(1,
                                                                           mean = theta.current[theta.estimated.names],
                                                                           sigma = covmat.proposal[theta.estimated.names,
                                                                                                   theta.estimated.names],
                                                                           lower = lower.proposal[theta.estimated.names],
                                                                           upper = upper.proposal[theta.estimated.names],
                                                                           algorithm = 'gibbs'))
    }
    target.theta.propose <- target(theta.propose)
    if (class(target.theta.propose) == "numeric") {
      suppressWarnings(target.theta.propose$log.density <- target.theta.propose)
      suppressWarnings(target.theta.propose$trace <- theta.propose)
    }
    if (!is.finite(target.theta.propose$log.density)) {
      log.acceptance <- -Inf
    }
    else {
      log.acceptance <- target.theta.propose$log.density -
        target.theta.current$log.density
      if(det(covmat.proposal)==0) {
        covmat.proposal <- covmat.proposal[theta.estimated.names,theta.estimated.names]+diag(ncol(covmat.proposal))*0.01 ### change made here
      }
      log.acceptance <- log.acceptance + tmvtnorm::dtmvnorm(x = theta.current[theta.estimated.names],
                                                            mean = theta.propose[theta.estimated.names],
                                                            sigma = covmat.proposal[theta.estimated.names,
                                                                                    theta.estimated.names], lower = lower.proposal[theta.estimated.names],
                                                            upper = upper.proposal[theta.estimated.names],
                                                            log = TRUE)
      log.acceptance <- log.acceptance - tmvtnorm::dtmvnorm(x = theta.propose[theta.estimated.names],
                                                            mean = theta.current[theta.estimated.names],
                                                            sigma = covmat.proposal[theta.estimated.names,
                                                                                    theta.estimated.names], lower = lower.proposal[theta.estimated.names],
                                                            upper = upper.proposal[theta.estimated.names],
                                                            log = TRUE)
    }
    if (verbose) {
      message(# "Propose: ", theta.propose[theta.estimated.names],
              "Proposal target: ", target.theta.propose[["log.density"]],
              ", acc prob: ", exp(log.acceptance), ", ", "acceptance rate: ",
              sprintf("%.3f", acceptance.rate), ", ", appendLF = FALSE)
    }
    if (is.accepted <- (log(runif(1)) < log.acceptance)) {
      theta.current <- theta.propose
      target.theta.current <- target.theta.propose
      if (verbose) {
        message("accepted")
      }
    }
    else if (verbose) {
      message("rejected")
    }
    trace <- rbind(trace, c(target.theta.current$trace))
    if (i.iteration == 1) {
      acceptance.rate <- is.accepted
    }
    else {
      if (is.null(acceptance.rate.weight)) {
        if (is.null(acceptance.window)) {
          acceptance.rate <- acceptance.rate + (is.accepted -
                                                  acceptance.rate)/i.iteration
        }
        else {
          acceptance.series <- c(is.accepted, acceptance.series)
          if (length(acceptance.series) > acceptance.window) {
            acceptance.series <- acceptance.series[-length(acceptance.series)]
          }
          acceptance.rate <- mean(acceptance.series)
        }
      }
      else {
        acceptance.rate <- acceptance.rate * (1 - acceptance.rate.weight) +
          is.accepted * acceptance.rate.weight
      }
    }
    tmp <- fitR::updateCovmat(covmat.empirical, theta.mean, theta.current,
                              i.iteration)
    covmat.empirical <- tmp$covmat
    theta.mean <- tmp$theta.mean

    if (check_variation('save_every_n') && i.iteration %% gc_testing$save_every_n == 0)
      saveRDS(trace, gc_env$partial_output_path)

  }

  return(list(trace = trace, acceptance.rate = acceptance.rate,
              covmat.empirical = covmat.empirical))
}
