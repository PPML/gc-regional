
new_ba <- readRDS("interventions_data_mobile_outreach_sensitivity_4_10th_BA.rds")
old_ba <- readRDS("interventions_data_mobile_outreach_sensitivity_BA.rds")

new_sf <- readRDS("interventions_data_mobile_outreach_sensitivity_4_10th_SF.rds")
old_sf <- readRDS("interventions_data_mobile_outreach_sensitivity_SF.rds")

library(magrittr)
library(dplyr)

old_ba %<>% filter(interv != '1year_high_activity_semi_annual_4_10th_hr_pt111_lr')
old_sf %<>% filter(interv != '1year_high_activity_semi_annual_4_10th_hr_pt111_lr')

ba <- rbind.data.frame(old_ba, new_ba)
sf <- rbind.data.frame(old_sf, new_sf)

saveRDS(ba, 'interventions_data_mobile_outreach_sensitivity_BA_fixed.rds')
saveRDS(sf, 'interventions_data_mobile_outreach_sensitivity_SF_fixed.rds')
