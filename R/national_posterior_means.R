#' Replace National Posterior Mean Values
#'

replace_national_parameter_medians <- function(theta) {
  if (! exists('national_posterior_medians', envir = gc_env)) {
    gc_env$national_posterior_medians <- apply(gc_env$national_natural_history, 2, median)
  }
  for (n in names(gc_env$national_posterior_medians)) {
    theta[[n]] <- gc_env$national_posterior_medians[[n]]
  }
  return(theta)
}

#' Remove Parameters which are Fixed for use with National Posterior Mean Values
remove_national_parameters <- function(theta) {
  theta <- theta[dplyr::setdiff(names(theta), names(gc_env$national_posterior_medians)) ]
  return(theta)
}
