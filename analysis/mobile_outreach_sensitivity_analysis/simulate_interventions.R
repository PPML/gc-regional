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

interventions_data <- simulate_interventions(
  interventions =
c(
'base_case',
'1year_high_activity_semi_annual',
'1year_high_activity_semi_annual_lower_capture_rate',
'1year_high_activity_semi_annual_complete_capture_rate',
'1year_high_activity_semi_annual_0_10th_hr_pt1555_lr',
'1year_high_activity_semi_annual_2_10th_hr_pt133_lr',
'1year_high_activity_semi_annual_4_10th_hr_pt111_lr',
'1year_high_activity_semi_annual_5_10th_hr_pt1_lr',
'1year_high_activity_semi_annual_6_10th_hr_pt088_lr',
'1year_high_activity_semi_annual_8_10th_hr_pt066_lr'), trace=trace)

saveRDS(interventions_data, paste0("interventions_data_mobile_outreach_sensitivity_", site, ".rds"))

print(paste0(site, " done!"))
