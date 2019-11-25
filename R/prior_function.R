###################################################
###   function to calculate prior likelihoods   ###
###################################################
#reads in current value for each parameter (contained in theta), and calculates likelihood
#asssumed distribution type for each parameter, and associated parameters describing each distribution contained in input file 'priors
#prior.param1 = 1st parameter for distribution; prior.param2=2nd param for distribution


prior_components2 <- function(theta) {
  stopifnot(length(theta) == nrow(gc_env$priors))

  vapply(1:length(theta), function(i) {
    gc_env$distribution_density_lookup[[gc_env$priors[[i, 'distribution']]]](
      gc_env$inverse_transforms[[strsplit(names(gc_env$theta)[[i]], "\\.")[[1]][[1]]]](theta[[i]]),
      gc_env$priors[[i, 'param1']],
      gc_env$priors[[i, 'param2']],
      log = T
    )
  },
  0
  )
}

# prior made from prior_components2
prior <- function(theta) sum(prior_components2(theta))


# Get the qth Quantile Value for Each of the Prior Distributions
qprior_transf <- function(q) {
  apply(gc_env$priors, 1, function(row) {
    gc_env$transforms[[row[['transformation']]]](
      gc_env$quantile_funs[[row[['distribution']]]](
        q,
        as.numeric(row[['param1']]),
        as.numeric(row[['param2']])
      )
    )
  })
}

rprior_transf <- function() {
  out <- apply(gc_env$priors, 1, function(row) {
    gc_env$transforms[[row[['transformation']]]](
      gc_env$random_funs[[row[['distribution']]]](
        1,
        as.numeric(row[['param1']]),
        as.numeric(row[['param2']])
      )
    )
  })
  names(out) <- names(gc_env$theta)
  return(out)
}



old_prior <- function(theta) {
  if (check_variation('imis')) {
    return(
      apply(theta, 1, function(x) exp(sum(prior_components(unlist(x)))))
      )
  } else {
    return(sum(prior_components(theta)))
  }
}

prior_components <- function(theta) {
  theta.names <- names(theta)
  priors <- gc_env$priors
  # Drop priors that aren't mentioned in theta
  priors <- dplyr::filter(priors, paste0(transformation, '.', parameter) %in% theta.names)

  # Construct a vector matching theta.names to the right rows of the priors df.
  # get prefixes before the first . in each theta.name
  prefixes <- sapply(strsplit(theta.names, "\\."), `[[`, 1)
  bad_prefixes <- ! prefixes %in% names(gc_env$inverse_transforms)

  # error if there are bad_prefixes
  if (any(bad_prefixes)) {
    x.names <- theta.names[which(bad_prefixes)]
    stop("The theta vector contains improperly formatted entries.\n",
         "The names '", paste(x.names, collapse = "', '"), "' start with prefixes \n",
         "which are not one of the transformations: \n",
         paste(as.character(names(gc_env$inverse_transforms)), collapse = ", "))
  }

  # get the parameters without the transformation prefix
  parameters <- gsub(
    pattern =
      paste(paste0(names(gc_env$inverse_transforms), "."), collapse="|"),
    replacement = "",
    x = theta.names
  )

  # for each parameter in theta, get the matching priors row
  parameters_lookup <- sapply(1:length(parameters), function(x) {
    which(priors$parameter == parameters[[x]])
  })

  # error if there are non-unique matches
  bad_matching_length <- sapply(parameters_lookup, length) != 1
  if (any(bad_matching_length)) {
    x.names <- theta.names[which(bad_matching_length)]
    stop("The following parts of theta do not uniquely match parameters listed\n",
         "in the priors dataframe:\n",
         paste(strwrap(paste(x.names, collapse = ", "), width = 80), collapse = "\n")
         )
  }

  # name the parameters_lookup vector
  names(parameters_lookup) <- theta.names

  # compute prior density in each component
  prior_vals <- sapply(1:length(theta), function(x) {

    #  ---- Error if x.name's prefix isn't a transformation ----
    # First let's confirm that the name of the theta component is formatted properly.
    # it should start with an inverse transformation.
    x.name <- theta.names[[x]]

    prior_density_f <-
      gc_env$distribution_density_lookup[[priors[[parameters_lookup[[x.name]], 'distribution']]]]

    return(
      prior_density_f(
        gc_env$inverse_transforms[[as.character(priors[x, 'transformation'])]](theta[[x]]),
        as.numeric(priors[x, 'param1']),
        as.numeric(priors[x, 'param2']),
        log = T)
    )
  })
  return(prior_vals)
}


prior_old <- function(theta) {
  prior.param1 <- gc_env$prior.param1
  prior.param2 <- gc_env$prior.param2
  national_natural_history <- gc_env$national_natural_history

  dbeta(ilogit(theta["logit.epsilon.1"]),prior.param1["epsilon.1"],prior.param2["epsilon.1"])*
  dbeta(ilogit(theta["logit.epsilon.2"]),prior.param1["epsilon.2"],prior.param2["epsilon.2"])*
  dbeta(ilogit(theta["logit.epsilon.3"]),prior.param1["epsilon.3"],prior.param2["epsilon.3"])*
  dbeta(ilogit(theta["logit.epsilon.4"]),prior.param1["epsilon.4"],prior.param2["epsilon.4"])*
  dgamma(exp(theta["log.rr.screen.m3"]),prior.param1["rr.screen.m3"],prior.param2["rr.screen.m3"]) *
  dgamma(exp(theta["log.rr.screen.f3"]),prior.param1["rr.screen.f3"],prior.param2["rr.screen.f3"]) *
  dgamma(exp(theta["log.rr.screen.ac"]),prior.param1["rr.screen.ac"],prior.param2["rr.screen.ac"]) *
  dgamma(exp(theta["log.c.min.m.1"]),prior.param1["c.min.m1"],prior.param2["c.min.m1"]) *
  dgamma(exp(theta["log.c.min.m.2"]),prior.param1["c.min.m2"],prior.param2["c.min.m2"]) *
  dgamma(exp(theta["log.c.min.f.1"]),prior.param1["c.min.f1"],prior.param2["c.min.f1"]) *
  dgamma(exp(theta["log.c.min.f.2"]),prior.param1["c.min.f2"],prior.param2["c.min.f2"])*
  dgamma(exp(theta["log.c.min.msm.1"]),prior.param1["c.min.msm1"],prior.param2["c.min.msm1"]) *
  dgamma(exp(theta["log.c.min.msm.2"]),prior.param1["c.min.msm2"],prior.param2["c.min.msm2"]) *
  dbeta(ilogit(theta["logit.b.m"]), prior.param1["b.m"],prior.param2["b.m"]) *
  dbeta(ilogit(theta["logit.b.f"]), prior.param1["b.f"],prior.param2["b.f"]) *
  dbeta(ilogit(theta["logit.b.msm"]), prior.param1["b.msm"],prior.param2["b.msm"])*
  infer_density_function(theta["log.dur.inf.symp.m"], national_natural_history$log.dur.inf.symp.m) *
  infer_density_function(theta["log.dur.inf.symp.msm"], national_natural_history$log.dur.inf.symp.msm) *
  infer_density_function(theta["log.dur.inf.symp.f"], national_natural_history$log.dur.inf.symp.f) *
  infer_density_function(theta["log.dur.inf.asymp.m"], national_natural_history$log.dur.inf.asymp.m) *
  infer_density_function(theta["log.dur.inf.asymp.f"], national_natural_history$log.dur.inf.asymp.f) *
  infer_density_function(theta["log.dur.inf.asymp.msm"], national_natural_history$log.dur.inf.asymp.msm) *
  infer_density_function(theta["logit.symp.m"], national_natural_history$logit.symp.m) *
  infer_density_function(theta["logit.symp.f"], national_natural_history$logit.symp.f) *
  infer_density_function(theta["logit.symp.msm"], national_natural_history$logit.symp.msm) *
  dbeta(ilogit(theta["logit.theta.1"]),prior.param1["theta.1"],prior.param2["theta.1"]) *
  dbeta(ilogit(theta["logit.theta.2"]),prior.param1["theta.2"],prior.param2["theta.2"]) *
  dbeta(ilogit(theta["logit.theta.3"]),prior.param1["theta.3"],prior.param2["theta.3"]) *
  dbeta(ilogit(theta["logit.theta.4"]),prior.param1["theta.4"],prior.param2["theta.4"])*
  dbeta(ilogit(theta["logit.theta.5"]),prior.param1["theta.5"],prior.param2["theta.5"]) *
  dbeta(ilogit(theta["logit.theta.6"]),prior.param1["theta.6"],prior.param2["theta.6"]) *
  dbeta(ilogit(theta["logit.theta.7"]),prior.param1["theta.7"],prior.param2["theta.7"]) *
  dbeta(ilogit(theta["logit.pi.m"]),prior.param1["pi.m"],prior.param2["pi.m"]) *
  dbeta(ilogit(theta["logit.pi.f"]),prior.param1["pi.f"],prior.param2["pi.f"]) *
  dbeta(ilogit(theta["logit.pi.msm"]),prior.param1["pi.msm"],prior.param2["pi.msm"]) *
  dnorm(exp(theta["log.rp.1.1.2.2"]),prior.param1["rp.1.1.2.2"], prior.param2["rp.1.1.2.2"])*
  dnorm(exp(theta["log.rp.4.1.2.1"]),prior.param1["rp.4.1.2.1"], prior.param2["rp.4.1.2.1"])*
  dnorm(exp(theta["log.rp.4.1.2.2"]),prior.param1["rp.4.1.2.2"], prior.param2["rp.4.1.2.2"])*
  dgamma(exp(theta["log.rp.1.1.1.1"]),prior.param1["rp.1.1.1.1"], prior.param2["rp.1.1.1.1"])*
  dnorm(exp(theta["log.rp.1.1.2.1"]),prior.param1["rp.1.1.2.1"], prior.param2["rp.1.1.2.1"])*
  dgamma(exp(theta["log.rp.1.2.1.1"]),prior.param1["rp.1.2.1.1"], prior.param2["rp.1.2.1.1"])*
  dgamma(exp(theta["log.rp.1.2.2.1"]),prior.param1["rp.1.2.2.1"], prior.param2["rp.1.2.2.1"])*
  dgamma(exp(theta["log.rp.1.1.1.2"]),prior.param1["rp.1.1.1.2"], prior.param2["rp.1.1.1.2"])*
  dgamma(exp(theta["log.rp.1.2.1.2"]),prior.param1["rp.1.2.1.2"], prior.param2["rp.1.2.1.2"])*
  dgamma(exp(theta["log.rp.1.2.2.2"]),prior.param1["rp.1.2.2.2"], prior.param2["rp.1.2.2.2"])*
  dgamma(exp(theta["log.rp.2.1.1.1"]),prior.param1["rp.2.1.1.1"], prior.param2["rp.2.1.1.1"])*
  dnorm(exp(theta["log.rp.2.1.2.1"]),prior.param1["rp.2.1.2.1"], prior.param2["rp.2.1.2.1"])*
  dgamma(exp(theta["log.rp.2.2.1.1"]),prior.param1["rp.2.2.1.1"], prior.param2["rp.2.2.1.1"])*
  dgamma(exp(theta["log.rp.2.2.2.1"]),prior.param1["rp.2.2.2.1"], prior.param2["rp.2.2.2.1"])*
  dgamma(exp(theta["log.rp.2.1.1.2"]),prior.param1["rp.2.1.1.2"], prior.param2["rp.2.1.1.2"])*
  dnorm(exp(theta["log.rp.2.1.2.2"]),prior.param1["rp.2.1.2.2"], prior.param2["rp.2.1.2.2"])*
  dgamma(exp(theta["log.rp.2.2.1.2"]),prior.param1["rp.2.2.1.2"], prior.param2["rp.2.2.1.2"])*
  dgamma(exp(theta["log.rp.2.2.2.2"]),prior.param1["rp.2.2.2.2"], prior.param2["rp.2.2.2.2"])*
  dgamma(exp(theta["log.rp.3.1.2.1"]),prior.param1["rp.3.1.2.1"], prior.param2["rp.3.1.2.1"])*
  dgamma(exp(theta["log.rp.3.2.2.1"]),prior.param1["rp.3.2.2.1"], prior.param2["rp.3.2.2.1"])*
  dgamma(exp(theta["log.rp.3.1.2.2"]),prior.param1["rp.3.1.2.2"], prior.param2["rp.3.1.2.2"])*
  dgamma(exp(theta["log.rp.3.2.2.2"]),prior.param1["rp.3.2.2.2"], prior.param2["rp.3.2.2.2"])*
  dbeta(ilogit(theta["logit.rep.symp.a"]),prior.param1["rep.symp.a"],prior.param2["rep.symp.a"]) *
  dbeta(ilogit(theta["logit.rep.symp.d"]),prior.param1["rep.symp.d"],prior.param2["rep.symp.d"])*
  dbeta(ilogit(theta["logit.rand.rep.symp.b"]),prior.param1["rand.rep.symp.b"],prior.param2["rand.rep.symp.b"])*
  dbeta(ilogit(theta["logit.rand.rep.symp.c"]),prior.param1["rand.rep.symp.c"],prior.param2["rand.rep.symp.c"])*
  dbeta(ilogit(theta["logit.screen.f1.a"]),prior.param1["screen.f1.a"],prior.param2["screen.f1.a"])*
  dbeta(ilogit(theta["logit.screen.f1.d"]),prior.param1["screen.f1.d"],prior.param2["screen.f1.d"])*
  dbeta(ilogit(theta["logit.rand.screen.f1.b"]),prior.param1["rand.screen.f1.b"],prior.param2["rand.screen.f1.b"])*
  dbeta(ilogit(theta["logit.rand.screen.f1.c"]),prior.param1["rand.screen.f1.c"],prior.param2["rand.screen.f1.c"])*
  dbeta(ilogit(theta["logit.screen.f2.a"]),prior.param1["screen.f2.a"],prior.param2["screen.f2.a"])*
  dbeta(ilogit(theta["logit.screen.f2.d"]),prior.param1["screen.f2.d"],prior.param2["screen.f2.d"])*
  dbeta(ilogit(theta["logit.rand.screen.f2.b"]),prior.param1["rand.screen.f2.b"],prior.param2["rand.screen.f2.b"])*
  dbeta(ilogit(theta["logit.rand.screen.f2.c"]),prior.param1["rand.screen.f2.c"],prior.param2["rand.screen.f2.c"])*
  dbeta(ilogit(theta["logit.screen.f1.1.a"]),prior.param1["screen.f1.1.a"],prior.param2["screen.f1.1.a"])*
  dbeta(ilogit(theta["logit.screen.f1.1.d"]),prior.param1["screen.f1.1.d"],prior.param2["screen.f1.1.d"])*
  dbeta(ilogit(theta["logit.rand.screen.f1.1.b"]),prior.param1["rand.screen.f1.1.b"],prior.param2["rand.screen.f1.1.b"])*
  dbeta(ilogit(theta["logit.rand.screen.f1.1.c"]),prior.param1["rand.screen.f1.1.c"],prior.param2["rand.screen.f1.1.c"])*
  dbeta(ilogit(theta["logit.screen.f2.1.a"]),prior.param1["screen.f2.1.a"],prior.param2["screen.f2.1.a"])*
  dbeta(ilogit(theta["logit.screen.f2.1.d"]),prior.param1["screen.f2.1.d"],prior.param2["screen.f2.1.d"])*
  dbeta(ilogit(theta["logit.rand.screen.f2.1.b"]),prior.param1["rand.screen.f2.1.b"],prior.param2["rand.screen.f2.1.b"])*
  dbeta(ilogit(theta["logit.rand.screen.f2.1.c"]),prior.param1["rand.screen.f2.1.c"],prior.param2["rand.screen.f2.1.c"])*
  dbeta(ilogit(theta["logit.screen.m1.a"]),prior.param1["screen.m1.a"],prior.param2["screen.m1.a"])*
  dbeta(ilogit(theta["logit.screen.m1.d"]),prior.param1["screen.m1.d"],prior.param2["screen.m1.d"])*
  dbeta(ilogit(theta["logit.rand.screen.m1.b"]),prior.param1["rand.screen.m1.b"],prior.param2["rand.screen.m1.b"])*
  dbeta(ilogit(theta["logit.rand.screen.m1.c"]),prior.param1["rand.screen.m1.c"],prior.param2["rand.screen.m1.c"])*
  dbeta(ilogit(theta["logit.screen.m2.a"]),prior.param1["screen.m2.a"],prior.param2["screen.m2.a"])*
  dbeta(ilogit(theta["logit.screen.m2.d"]),prior.param1["screen.m2.d"],prior.param2["screen.m2.d"])*
  dbeta(ilogit(theta["logit.rand.screen.m2.b"]),prior.param1["rand.screen.m2.b"],prior.param2["rand.screen.m2.b"])*
  dbeta(ilogit(theta["logit.rand.screen.m2.c"]),prior.param1["rand.screen.m2.c"],prior.param2["rand.screen.m2.c"])*
  dbeta(ilogit(theta["logit.screen.m1.1.a"]),prior.param1["screen.m1.1.a"],prior.param2["screen.m1.1.a"])*
  dbeta(ilogit(theta["logit.screen.m1.1.d"]),prior.param1["screen.m1.1.d"],prior.param2["screen.m1.1.d"])*
  dbeta(ilogit(theta["logit.rand.screen.m1.1.b"]),prior.param1["rand.screen.m1.1.b"],prior.param2["rand.screen.m1.1.b"])*
  dbeta(ilogit(theta["logit.rand.screen.m1.1.c"]),prior.param1["rand.screen.m1.1.c"],prior.param2["rand.screen.m1.1.c"])*
  dbeta(ilogit(theta["logit.screen.m2.1.a"]),prior.param1["screen.m2.1.a"],prior.param2["screen.m2.1.a"])*
  dbeta(ilogit(theta["logit.screen.m2.1.d"]),prior.param1["screen.m2.1.d"],prior.param2["screen.m2.1.d"])*
  dbeta(ilogit(theta["logit.rand.screen.m2.1.b"]),prior.param1["rand.screen.m2.1.b"],prior.param2["rand.screen.m2.1.b"])*
  dbeta(ilogit(theta["logit.rand.screen.m2.1.c"]),prior.param1["rand.screen.m2.1.c"],prior.param2["rand.screen.m2.1.c"])*
  dbeta(ilogit(theta["logit.screen.msm1.a"]),prior.param1["screen.msm1.a"],prior.param2["screen.msm1.a"])*
  dbeta(ilogit(theta["logit.screen.msm1.d"]),prior.param1["screen.msm1.d"],prior.param2["screen.msm1.d"])*
  dbeta(ilogit(theta["logit.rand.screen.msm1.b"]),prior.param1["rand.screen.msm1.b"],prior.param2["rand.screen.msm1.b"])*
  dbeta(ilogit(theta["logit.rand.screen.msm1.c"]),prior.param1["rand.screen.msm1.c"],prior.param2["rand.screen.msm1.c"])*
  dbeta(ilogit(theta["logit.screen.msm2.a"]),prior.param1["screen.msm2.a"],prior.param2["screen.msm2.a"])*
  dbeta(ilogit(theta["logit.screen.msm2.d"]),prior.param1["screen.msm2.d"],prior.param2["screen.msm2.d"])*
  dbeta(ilogit(theta["logit.rand.screen.msm2.b"]),prior.param1["rand.screen.msm2.b"],prior.param2["rand.screen.msm2.b"])*
  dbeta(ilogit(theta["logit.rand.screen.msm2.c"]),prior.param1["rand.screen.msm2.c"],prior.param2["rand.screen.msm2.c"])*
  dbeta(ilogit(theta["logit.behav.lin"]),prior.param1["behav.lin"],prior.param2["behav.lin"])*
  dbeta(ilogit(theta["logit.risk.rep.symp.m"]),prior.param1["risk.rep.symp.m"],prior.param2["risk.rep.symp.m"]) *
  dbeta(ilogit(theta["logit.risk.rep.symp.m1"]),prior.param1["risk.rep.symp.m1"],prior.param2["risk.rep.symp.m1"]) *
  dbeta(ilogit(theta["logit.risk.rep.symp.f"]),prior.param1["risk.rep.symp.f"],prior.param2["risk.rep.symp.f"]) *
  dbeta(ilogit(theta["logit.risk.rep.symp.f1"]),prior.param1["risk.rep.symp.f1"],prior.param2["risk.rep.symp.f1"]) *
  ( if (logistic_screening_variation()) {
    prod(
      dgamma(
      x = sapply(unique(gc_env$screening_names), function(x) {
      exp(theta[paste0('log.', x, '.growth')])
      })
      ,shape = 2
      ,rate = 3
      )
    ) *
    dgamma(exp(theta['log.rep.symp.growth']), 2, 3)
  } else 1
  )
}


threshold_distribution <- function(x, min, max, value) {
  if (x<min) { return(0) }
  if (x>max) { return(0) }
  return(value)
}


infer_density_function <- function(x, vec) {
  stopifnot(is.numeric(vec))
  d <- density(vec)

  # Add a top and bottom to the density values which are 0, so that any value
  # outside the support of the density vector is 0 probability.
  step <- abs(d$x[[2]] - d$x[[1]])
  bottom_bin <- d$x[[1]] - step
  d$x <- append(x = d$x, values = bottom_bin, after = 0)
  d$y <- append(x = d$y, values = 0, after = 0)
  top_bin <- d$x[[length(d$x)]] + step
  d$x <- append(x = d$x, values = top_bin, after = length(d$x))
  d$y <- append(x = d$y, values = 0, after = length(d$y))

  # lookup closest neighbor to x
  dist_x <- abs(d$x - x)
  closest_x_index <- which(dist_x == min(dist_x))[[1]]

  # get y value
  return(d$y[[closest_x_index]])
}

