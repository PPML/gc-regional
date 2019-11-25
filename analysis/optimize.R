# Optimize gcRegional
#
#    In this script we're going to set the site, specify the variation of the
#    model, and to sample the priors randomly to start off several chains of
#    calibrations.
#
#    - Set up gcRegional
#    - Write qprior_transf
#    - Write rprior_transf
#    - Write transform_theta
#    - Rewrite the Prior
#    - Sample the Prior
#    - Run Optimization Loop
#
#    This script will be used on the cluster, run 50 times in parallel,
#    (25 for each site), and will start from a random place, so that we can
#    see if the optimizations converge the same place.
#

# Script Config:
library(magrittr)
library(gcRegional)

N_loops <- 40


site <-
  if (on_the_cluster()) {
    commandArgs(trailingOnly = T)[[1]]
  } else "SF" # default to one as a test
# and use commandLineArgs() if we're running
# on the cluster.

# minimum dLogPosterior acceptable for starting an optimization chain:
min_acceptable <- -10000
nth_sim <- if (on_the_cluster()) {
  nth_sim <-as.numeric(commandArgs(trailingOnly=T)[[2]])
} else 1
set.seed(nth_sim)

output_directory <-
  if(on_the_cluster()) "~/2019/January/22-2/optim/" else "~/Documents/gcRegional/output/09-23-18/optim/"

# Set up gcRegional -------------------------------------------------------
#
# - Specify our model variation:
#
#   We're weighting our beta distributions with denominators scaled by sample
#   sizes in the SSuN case report data.  We're going to try scaling the NHANES
#   population size down by a factor of 100 and the SSuN Case reports by 10.
#
# - Configure the model

gc_testing <- c(
  'population_denominators',
  'msm_incidence_likelihood',
  'increase_risk_all_pops',
  'separate_increased_risk',
  'only_use_nhanes_women'
)

rm(gc_env) # remove any data lurking from a previous simulation

load_start(site)

gc_env$scalars <- c(nhanes_prev = 1/100, case_rates = 1/10)


theta_list <- list()
optim_list <- list()


# Run an Optimization Loop ------------------------------------------------
#
# - Randomly Sample the Priors
# - Loop:
#     - Run Optim on dLogPosterior with Sample
#     - Store Output after Each Optim
#     - Set up the next to start with the output of optim
#

# Randomly Sample the Priors until Acceptable
sample_until_acceptable <- function() {
  theta <- rprior_transf()
  accepted = F
  while(! accepted) {
    theta <- rprior_transf()
    d <- dLogPosterior(theta)
    print(d)
    accepted = is.finite(d) && d > min_acceptable
  }
  return(theta)
}
theta <- sample_until_acceptable()

# Optimization Loop

if (! on_the_cluster()) {
  dLogPosterior <- function(x) {
    x <- gcRegional::dLogPosterior(x)
    print(x)
    return(x)
  }
}

n_done <- 0
done <- F
while (!done) {
  out <- optim(
    par = theta,
    fn = dLogPosterior,
    method = c("BFGS", "Nelder-Mead")[[if (n_done > 10) 1 else 2]],
    control = list(fnscale = -1, maxit=1000)
  )

  if (out$value != -1e32) {
    n_done <- n_done + 1
    done = n_done == N_loops

    theta <- out$par
    names(theta) <- names(gc_env$theta)
    theta_list[[n_done]] <- theta
    optim_list[[n_done]] <- out
    saveRDS(optim_list,
            paste0(
              output_directory,
              site,
              "_optim_list_",
              nth_sim,
              ".rds"
            )
    )
  } else {
    theta <- sample_until_acceptable()
  }
}
