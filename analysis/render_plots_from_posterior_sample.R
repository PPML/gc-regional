devtools::load_all(".")
library(coda)
library(fitR)

gc_testing <- c(
  'population_denominators', 'msm_incidence_likelihood', 'increase_risk_all_pops', 'separate_increased_risk', 'only_use_nhanes_women'
)


for (site in c('SF', 'BA')) {
  load_start(site)
  output_directory <-
    if (on_the_cluster()) "~/2018/December/17/mcmc/" else "~/Documents/gcRegional/output/12-13-18/mcmc/"
  post.sample <- readRDS(paste0(output_directory, site, "_posterior_sample.rds"))
  trace.burn.thin <- trace.burn <- trace <- as.data.frame(post.sample$theta.list)
  pred <- as.data.frame(post.sample$outputs)
  gc_assign(pred)
  plot_posteriors(output_directory)
}

