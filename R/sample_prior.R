distribution_random_lookup <- list(
  beta = rbeta,
  gamma = rgamma,
  normal = rnorm
)

transformation_lookup <- list(
  logit = logit,
  log = log
)

sample_prior <- function() {
  priors <- gc_env$priors
  sample <- suppressWarnings(apply(priors, 1, function(x) {
    transformation_lookup[[x[['transformation']]]](
      distribution_random_lookup[[x[['distribution']]]](1, as.numeric(x[['param1']]), as.numeric(x[['param2']]))
    )
  }))
  nan_entries <- which(is.nan(sample))
  # Since some of the priors are log transformed but use normal distributions
  # there is positive probability of computing the log of a negative number
  # in the above computation. To correct for this, we will find the NaN
  # entries and sample the priors again until they are acceptable.
  while (length(nan_entries) > 0) {
    for (nan_entry in nan_entries) {
      x <- priors[nan_entry, ]
      sample[[nan_entry]] <-
        transformation_lookup[[x[['transformation']]]](
          distribution_random_lookup[[x[['distribution']]]](
            1, as.numeric(x[['param1']]), as.numeric(x[['param2']])
          )
        )
    }
    nan_entries <- which(is.nan(sample))
  }
  names(sample) <- paste0(priors$transformation, ".", priors$parameter)
  return(sample)
}
