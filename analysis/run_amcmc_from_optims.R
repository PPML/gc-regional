# Run AMCMC With Top Values from Optimization
library(adaptMCMC)
library(gcRegional)

if (on_the_cluster()) {
  args <- commandArgs(trailingOnly = T)
  site <- args[[1]]
  array_id <- args[[2]]
  gcR_path <- "~/gcRegional/"
  output_dir <- "~/2018/December/13/mcmc/"

} else {
  site <- c("SF", "BA")[[rbinom(1, 1, .5)+1]]
  array_id <- 1
  gcR_path <- "~/Documents/gcRegional/"
  output_dir <- "~/Documents/gcRegional/output/11-05-18/"
}

gc_testing <-
  c(
    'population_denominators',
    'msm_incidence_likelihood',
    'increase_risk_all_pops',
    'separate_increased_risk',
    'only_use_nhanes_women',
    'only_increase_hetero_high_activity_transmission'
  )
load_start(site)
gc_env$scalars <- c(nhanes_prev = 1/100, case_rates = 1/10)

top_optim_trace_path <- paste0(gcR_path, "output/01-07-19/", site, "_optim_trace_top.rds")
optim_trace <- readRDS(top_optim_trace_path)

theta <- unlist(optim_trace[ (as.numeric(array_id) %% 15)+1, ])

sd_optims <- apply(optim_trace, 2, sd)

sample <- MCMC(
  p = dLogPosterior,
  n = 25000,
  init = theta,
  scale = sd_optims / 1000,
  adapt = T,
  acc.rate = .234, showProgressBar = T)

saveRDS(sample, file = paste0(output_dir, site, "_mcmc_25000_", array_id, ".rds"))
