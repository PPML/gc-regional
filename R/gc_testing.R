
logistic_screening_variation <- function() {
  if (exists('gc_testing', envir=.GlobalEnv) &&
      'logistic_screening' %in% gc_testing) {
    TRUE
  } else FALSE
}


only_increasing_screening <- function() {
  if (exists('gc_testing', envir=.GlobalEnv) &&
      'only_increasing_screening' %in% gc_testing) {
    TRUE
  } else FALSE
}


lower_proportion_msm <- function() {
  if (exists('gc_testing', envir=.GlobalEnv) &&
      'lower_proportion_msm' %in% gc_testing) {
    TRUE
  } else FALSE
}

widen_likelihood <- function() {
  if (exists("gc_testing", envir = .GlobalEnv) &&
      'widen_likelihood' %in% gc_testing) {
    TRUE
  } else FALSE
}

check_variation <- function(str) {
  if (exists('gc_testing', envir = .GlobalEnv) &&
      str %in% gc_testing) {
    TRUE
  } else FALSE
}
