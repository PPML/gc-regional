library(IMIS)

safeLikelihood <- function(theta) {
  tryCatch({
    ll <- gcRegional::likelihood(theta)$LogLL
    print(ll)
    return(ll)
  }, error = function() -1e32)
}

safePrior <- function(theta) {
  tryCatch({
    lp <- gcRegional::prior(theta)
    print(lp)
    return(lp)
  }, error = function() -1e32)
}

likelihood <- function(trace) {
  apply(trace, 1, function(theta) exp(safeLikelihood(unlist(theta))))
}

prior <- function(trace) {
  apply(trace, 1, function(theta) exp(safePrior(unlist(theta))))
}

sample.prior <- function(n) {
  return(optim_trace)
}

myIMIS(500, 3000, 100, 10)
