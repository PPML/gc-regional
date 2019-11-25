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

N_loops <- 20


site <-
  if (on_the_cluster()) {
    commandArgs(trailingOnly = T)[[1]]
  } else "SF" # default to one as a test
# and use commandLineArgs() if we're running
# on the cluster.

# minimum dLogPosterior acceptable for starting an optimization chain:
min_acceptable <- -4000
nth_sim <- if (on_the_cluster()) {
  nth_sim <-as.numeric(commandArgs(trailingOnly=T)[[2]])
} else 1
set.seed(nth_sim)

output_directory <-
  if(on_the_cluster()) "~/2018/December/5/optim/" else "~/Documents/gcRegional/output/11-02-18/optim/"

previous_optim_directory <-
  if(on_the_cluster()) "~/2018/November/28/optim/" else "~/Documents/gcRegional/output/11-01-18/optim-w-medians/"


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

gc_env$national_posterior_medians <- apply(gc_env$national_natural_history, 2, median)

theta_list <- list()
optim_list <- list()


# Randomly Sample the Priors until Acceptable
sample_until_acceptable <- function() {
  theta <- rprior_transf()
  accepted = F
  while(! accepted) {
    theta <- rprior_transf()
    d <- dLogPosterior(theta)
    accepted = is.finite(d) && d > min_acceptable
  }
  return(theta)
}

previous_optim <- readRDS(file.path(previous_optim_directory, paste0(site, "_optim_list_", nth_sim, ".rds")))
# theta <- previous_optim[[length(previous_optim)]][['par']]

previous_vals <- sapply(1:length(previous_optim), function(i) previous_optim[[i]][['value']])
max_val <- which(previous_vals == max(previous_vals))
theta <- previous_optim[[max_val[[1]]]][['par']]

if (check_variation('population_risk_mixing') && ! 'logit.epsilon' %in% names(theta))
  theta['logit.epsilon'] <- logit(rbeta(1,1.1,1.1))

if (check_variation('separate_increased_risk') && ! 'logit.behav.hetero' %in% names(theta))
  theta['logit.behav.hetero'] <- theta['logit.behav.lin']

if (! check_variation('national_posterior_medians') &&
    ! all(names(gc_env$national_posterior_medians) %in% names(theta))) {
  theta <- replace_national_parameter_medians(theta)
}

parameters <-  paste(gc_env$priors$transformation, gc_env$priors$parameter, sep =  '.')
theta <- sapply(parameters, function(x) theta[[x]])


# Optimization Loop

n_done <- 0
done <- F
while (!done) {
  out <- optim(
    par = theta,
    fn = dLogPosterior,
    method = c("BFGS", "Nelder-Mead")[[(n_done %% 2) + 1]],
    control = list(fnscale = -1, maxit=500)
  )

  if (out$value != -1e32) {
    n_done <- n_done + 1
    done = n_done == N_loops

    if (out$value > dLogPosterior(theta)) {
      theta <- out$par
      names(theta) <- names(gc_env$theta)
    }
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

