args <- commandArgs(trailingOnly=T)

site <- as.character(args[[1]])
print(site)

devtools::load_all()
library(magrittr)

gc_testing <- c(
'population_denominators',
'msm_incidence_likelihood',
'increase_risk_all_pops',
'separate_increased_risk',
'only_use_nhanes_women')
load_start(site)

trace <- switch(
	site,
	BA = as.data.frame(readRDS(system.file("calibration_outcomes/BA_posterior_sample.rds", package = "gcRegional"))$theta.list),
	SF = as.data.frame(readRDS(system.file("calibration_outcomes/SF_posterior_sample.rds", package = "gcRegional"))$theta.list)
)

interventions_data <- simulate_interventions()

saveRDS(interventions_data, paste0("interventions_data_", site, ".rds"))

print(paste0(site, " done!"))
