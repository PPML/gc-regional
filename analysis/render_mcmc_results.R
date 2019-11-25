# library(gcRegional)
devtools::load_all(".")
library(coda)
library(fitR)

gc_testing <- c(
  'population_denominators',
  'msm_incidence_likelihood',
  'increase_risk_all_pops',
  'separate_increased_risk',
  'only_use_nhanes_women',
  'only_increase_hetero_high_activity_transmission'
)

output_directory <-
  if (on_the_cluster()) "~/2018/November/29/mcmc/" else "~/Documents/gcRegional/output/01-22-19/mcmc/"

for (site in c("SF", "BA")) {
  rm(gc_env)
  load_start(site)
  calibration_files <- grep(paste0(site, "_mcmc"), list.files(output_directory, full.names = T), value = T)
  mcmc_list <- mcmc.list(lapply(calibration_files, function(f) mcmc(readRDS(f)$samples)))
  mcmc_list <- burnAndThin(mcmc_list, burn = 15000, thin = 100)
  trace <- do.call(rbind.data.frame, mcmc_list)
  print(paste0("nrow of trace.burn.thin = ", nrow(trace)))

  trace <- dplyr::sample_n(trace, 1000)

  trace.burn.thin <- trace.burn <- trace
  post.sample <- model_fits(trace, use_trace_without_sampling = T)
  pred <- as.data.frame(post.sample$outputs)
  gc_assign(pred)
  saveRDS(post.sample, paste0(output_directory, site, "_posterior_sample.rds"))
  plot_posteriors(output_directory)
}
