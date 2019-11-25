#' Combine, Clean, and Save Raw Interventions Simulation Data into interventions.rds 
#'
.render_interventions_df <- function() { 
  require(tidyverse)
  require(ggrepel)
  require(reshape2)
  require(magrittr)


  sf_df <- readRDS(
    system.file("intervention_analysis/intervention_data_SF3.rds", 
    package = 'gcRegional')) 
  ba_df <- readRDS(
  system.file("intervention_analysis/intervention_data_BA3.rds", 
  package = 'gcRegional')) 

  df <- do.call(rbind.data.frame, list(
    cbind.data.frame(site = 'SF', sf_df),
    cbind.data.frame(site = 'BA', ba_df)
  ))
  rm(sf_df)
  rm(ba_df)

  # Set the values from the 1st week of the year to be exactly 
  # yearly values. I.e. round down from dates like 2009.019 to 2009. 
  # Note that .019 ~ 1/52.
  # df$year <- sapply(df$year, function(x) if (x %% 1 < .25) x - x%%1 else x)
  df$year <- ifelse(df$year %% 1 < .25, df$year - df$year%%1, df$year)
  df$site <- as.character(df$site)
  df$interv <- as.character(df$interv)
  df$population <- as.character(df$population)


  # We didn't end up looking at the quarter annual interventions 
  df %<>% filter(! interv %in% c('no_screening', 'universal_annual',
  															 'universal_semi_annual'
  															 ))

  saveRDS(df, file.path(
    system.file("intervention_analysis/", package='gcRegional'),
    "interventions3.rds"))
}

#' Analyze interventions.df to Render and Save a Comparison DataFrame for Incidence vs. Diagnoses
#'
#' The data computed and saved here is used for plot rendered by 
#' plot_comparison_of_model_incidence_and_diagnoses() which appears 
#' in our supplemental materials.
.render_incidence_vs_reported_cases_data <- function() { 
  df <- readRDS(system.file("intervention_analysis/interventions3.rds", package='gcRegional'))

  last_data_year_sf <- 2016
  last_data_year_ba <- 2017


  # Subset the data frame for just incident_cases in the last data year and one year
  # prior, only populations of interest for the next plot, and just get the base_case.
  df %>% 
    select(site, sim, interv, year, population, incident_cases) %>% 
    filter(
      interv == 'base_case',
      population %in% c('Black Male', 'Black Female', 'Other Male', 'Other Female', 'Hispanic Male', 'Hispanic Female', 'MSM'),
      (year %in% c(last_data_year_sf, last_data_year_sf-1) & site == 'SF') |
        (year %in% c(last_data_year_ba, last_data_year_ba-1) & site == 'BA')
    ) %>% 
    select(-interv) -> cases_df


  # The incident_cases column is cumulative, so we will subtract the previous year from 
  # the next year. We split up the operation for SF and BA because the last data years are
  # different for each site.
  cases_df <- do.call(rbind.data.frame, list(
    dcast(filter(cases_df, site == 'SF'), site + sim + population ~ year, value.var = 'incident_cases') %>% 
      mutate(incidents = `2016` - `2015`) %>% 
      select(-c(`2016`, `2015`)),
    dcast(filter(cases_df, site == 'BA'), site + sim + population ~ year, value.var = 'incident_cases') %>% 
      mutate(incidents = (`2017` - `2016`)) %>% 
      select(-c(`2017`, `2016`))
  ))

  cases_df %>% 
    group_by(site, sim) %>% 
    mutate(fraction_of_all_incidents = incidents / sum(incidents)) -> fractions_of_incident_cases

  summary_cases_df <- group_by(fractions_of_incident_cases, site, population) %>% 
    summarize(
      ci_high = quantile(fraction_of_all_incidents, 0.975,  na.rm = T),
      iq_high = quantile(fraction_of_all_incidents, 0.75,  na.rm = T),
      mean = mean(fraction_of_all_incidents, na.rm = T),
      iq_low = quantile(fraction_of_all_incidents, 0.25,  na.rm = T),
      ci_low = quantile(fraction_of_all_incidents, 0.025,  na.rm = T))

  summary_cases_df$site <- sapply(as.character(summary_cases_df$site), function(x) c(BA = 'Baltimore', SF = 'San Francisco')[[x]])

  write.csv(x = summary_cases_df, file = file.path(
    system.file(
      package = 'gcRegional',
      "intervention_analysis/data/"
    ),
    'breakdown_of_case_distribution3.csv'
  ))


  sf_p_msm <- 
      read.delim(system.file(
        "extdata",
        paste("ssun_p_msm_", "SF", ".txt", sep = ""),
        package = 'gcRegional'
      ))


  ba_p_msm <- 
      read.delim(system.file(
        "extdata",
        paste("ssun_p_msm_", "BA", ".txt", sep = ""),
        package = 'gcRegional'
      ))



  sf_cases <- 
  tibble::tribble(
        ~sex,     ~age,      ~race, ~cases,
        "F", "15-24",    "black",       62,
        "F", "25-39",    "black",       67,
        "M", "15-24",    "black",       93,
        "M", "25-39",    "black",      221,
        "F", "15-24", "hispanic",       27,
        "F", "25-39", "hispanic",       30,
        "M", "15-24", "hispanic",      127,
        "M", "25-39", "hispanic",      462,
        "F", "15-24",    "other",       78,
        "F", "25-39",    "other",      158,
        "M", "15-24",    "other",      319,
        "M", "25-39",    "other",     1971
    )


  ba_cases <- 
  tibble::tribble(
    ~sex,   ~age,      ~race, ~cases, 
    "M", "15-24",    "black",    608, 
    "M", "25-39",    "black",    723, 
    "M", "15-24", "hispanic",      8, 
    "M", "25-39", "hispanic",     18, 
    "M", "15-24",  	 "other",    207, 
    "M", "25-39",	 	 "other",    219, 
    "F", "15-24",    "black",    831, 
    "F", "25-39",    "black",    390, 
    "F", "15-24", "hispanic",      6, 
    "F", "25-39", "hispanic",      4, 
    "F", "15-24",    "other",    293, 
    "F", "25-39",    "other",    168
    )

  reported_cases <- rbind.data.frame(
    cbind.data.frame(site = 'SF', sf_cases),
    cbind.data.frame(site = 'BA', ba_cases))

  sf_p_msm_y <- sf_p_msm %>% filter(year == 2016, age_cat==1) %>% `[[`('p_msm')
  sf_p_msm_o <- sf_p_msm %>% filter(year == 2016, age_cat==2) %>% `[[`('p_msm')
  ba_p_msm_y <- ba_p_msm %>% filter(year == 2017, age_cat==1) %>% `[[`('p_msm')
  ba_p_msm_o <- ba_p_msm %>% filter(year == 2017, age_cat==2) %>% `[[`('p_msm')

  p_msm_site_age <- tibble::tribble(
    ~site, ~sex, ~age, ~p_msm,
    'SF', 'F', '15-24', 0,
    'SF', 'F', '25-39', 0,
    'SF', 'M', '15-24', sf_p_msm_y,
    'SF', 'M', '25-39', sf_p_msm_o,
    'BA', 'F', '15-24', 0,
    'BA', 'F', '25-39', 0,
    'BA', 'M', '15-24', ba_p_msm_y,
    'BA', 'M', '25-39', ba_p_msm_o
  )

  reported_cases %<>% merge(p_msm_site_age)

  reported_cases %<>% mutate(
    cases_in_msw = cases * (1-p_msm),
    cases_in_msm = cases * p_msm
    )

  reported_cases %<>% select(site, sex, age, race, cases_in_msw, cases_in_msm) %>% 
    group_by(site, sex, race) %>% 
    summarize_if(is.numeric, sum)

  reported_cases %<>% mutate(population = paste(str_to_title(race), c('F' = 'Female', 'M' = 'Male')[sex])) %>% 
    select(-c(sex, race))

  reported_cases %<>% ungroup() %>% select(-sex)

  reported_cases_msw <- reported_cases[,c('site', 'population', 'cases_in_msw')]
  colnames(reported_cases_msw)[3] <- 'cases'

  reported_cases_msm <- tibble::tribble(
    ~site, ~population, ~cases,
    'SF', 'MSM', reported_cases %>% filter(site == 'SF') %>% `[[`('cases_in_msm') %>% sum(),
    'BA', 'MSM', reported_cases %>% filter(site == 'BA') %>% `[[`('cases_in_msm') %>% sum()
    )

  reported_cases <- rbind.data.frame(reported_cases_msw, reported_cases_msm)

  reported_cases %<>% mutate(site = c('SF' = 'San Francisco', 'BA' = 'Baltimore')[site])

  reported_cases %<>% group_by(site) %>% mutate(proportion_of_cases = cases / sum(cases))

  reported_cases %<>% select(-cases) %>% mutate(mean = proportion_of_cases,
                                                ci_high = mean, ci_low = mean,
                                                iq_high = mean, iq_low=mean) %>% select(-proportion_of_cases)

  # compare_reported_cases_and_estimated_incidents_proportions <- 
  diag_inc_comparison_df <- 
    do.call(rbind.data.frame, list(
      cbind.data.frame(variable = 'model_estimate', summary_cases_df), 
      cbind.data.frame(variable = 'reported_cases', reported_cases)))

  diag_inc_comparison_df_errorbars <- diag_inc_comparison_df
  diag_inc_comparison_df_errorbars[diag_inc_comparison_df_errorbars$variable == 'reported_cases',
    c('ci_high', 'iq_high', 'mean', 'iq_low', 'ci_low')] <- NA 


  saveRDS(diag_inc_comparison_df, file.path(
    system.file("intervention_analysis/data/", package='gcRegional'),
    "diag_inc_comparison_df3.rds"))

  saveRDS(diag_inc_comparison_df_errorbars, file.path(
    system.file("intervention_analysis/data/", package='gcRegional'),
    "diag_inc_comparison_df_errorbars3.rds"))
}


#' Construct interventions_descriptions so it is available in the Global Scope
.construct_interventions_description_df <- function() { 

  # df <- readRDS(system.file("intervention_analysis/interventions.rds", package='gcRegional'))

  pal1 <- RColorBrewer::brewer.pal(11, 'BrBG')
  pal2 <- RColorBrewer::brewer.pal(11, 'RdBu')
  pal3 <- RColorBrewer::brewer.pal(11, 'PiYG')
  pal4 <- sapply(1:11, function(x) c(pal1[x], pal2[x], pal3[x]))
  pal4 <- c(pal4[,-c(5,6,7)])[1:23]

  # > pal4
  #  [1] "#543005" "#67001F" "#8E0152" "#8C510A" "#B2182B" "#C51B7D" "#BF812D" "#D6604D"
  #  [9] "#DE77AE" "#DFC27D" "#F4A582" "#F1B6DA" "#80CDC1" "#92C5DE" "#B8E186" "#35978F"
  # [17] "#4393C3" "#7FBC41" "#01665E" 
  # ltfu + van: "#2166AC" "#4D9221" "#003C30" "#053061"

  intervention_descriptions <<- tibble::tribble(
                                    ~shortnames,         ~names,       ~frequency,    ~color,
                                    "base_case",         "Base Case",      "base_case", "dimgrey",

                                 "young_annual",         "15-24 Annual Screening",         "annual", "#543005",
                            "young_semi_annual",         "15-24 Twice-Annual Screening",    "semi_annual", "#8E0152",

                          "female_young_annual",         "Female 15-24 Annual Screening",         "annual", "#BF812D",
                     "female_young_semi_annual",         "Female 15-24 Twice-Annual Screening",    "semi_annual", "#D6604D",

                                   "MSM_annual",         "MSM Annual Screening",         "annual", "#B8E186",
                              "MSM_semi_annual",         "MSM Twice-Annual Screening",    "semi_annual", "#35978F",
                           "MSM_quarter_annual",         "MSM Quarter-Annual Screening", "quarter_annual", "#4393C3",

                   "high_activity_young_annual",         "High Risk 15-24 Annual Screening",         "annual", "#543005",
              "high_activity_young_semi_annual",         "High Risk 15-24 Twice-Annual Screening",    "semi_annual", "#67001F",

            "high_activity_female_young_annual",         "High Risk 15-24 Female Annual Screening",         "annual", "#BF812D",
       "high_activity_female_young_semi_annual",         "High Risk 15-24 Female Twice-Annual Screening",    "semi_annual", "#D6604D",

                     "high_activity_MSM_annual",         "High Risk MSM Annual Screening",         "annual", "#B8E186",
                "high_activity_MSM_semi_annual",         "High Risk MSM Twice-Annual Screening",    "semi_annual", "#35978F",
             "high_activity_MSM_quarter_annual",         "High Risk MSM Quarter-Annual Screening", "quarter_annual", "#4393C3",

                         "high_activity_annual",         "High Risk Annual Screening",         "annual", "#7FBC41",
                    "high_activity_semi_annual",         "High Risk Twice-Annual Screening",    "semi_annual", "#01665E",

                            "remove_ltfu_10pct",         "Remove 10% LTFU",           "ltfu", "#2166AC",
                            "remove_ltfu_20pct",         "Remove 20% LTFU",           "ltfu", "#4D9221",
        "remove_f_msm_10pct_and_msw_20pct_ltfu",         "Remove 10% LTFU (20% MSW)",           "ltfu", "#003C30",
              "1year_high_activity_semi_annual",         "Mobile Outreach Testing",    "semi_annual", "#053061",
"1year_high_activity_semi_annual_lower_capture_rate",    "Mobile Outreach Testing (25% HR, 12.78% LR)",    "semi_annual", "#053061",
'1year_high_activity_semi_annual_complete_capture_rate', "Mobile Outreach Testing (100% HR, .0444% LR)",    "semi_annual", "#053061",
'1year_high_activity_semi_annual_0_10th_hr_pt1555_lr',   "Mobile Outreach Testing (0% HR, 15.56% LR)",    "semi_annual", "#053061",
'1year_high_activity_semi_annual_2_10th_hr_pt133_lr',    "Mobile Outreach Testing (20% HR, 13.33% LR)",    "semi_annual", "#053061",
'1year_high_activity_semi_annual_4_10th_hr_pt111_lr',    "Mobile Outreach Testing (40% HR, 11.11% LR)",    "semi_annual", "#053061",
'1year_high_activity_semi_annual_5_10th_hr_pt1_lr',      "Mobile Outreach Testing (50% HR, 10.0% LR)",    "semi_annual", "#053061",
'1year_high_activity_semi_annual_6_10th_hr_pt088_lr',    "Mobile Outreach Testing (60% HR, 8.89% LR)",    "semi_annual", "#053061",
'1year_high_activity_semi_annual_8_10th_hr_pt066_lr',    "Mobile Outreach Testing (80% HR, 6.67% LR)",    "semi_annual", "#053061"
    )


  empirical_interventions <<- c('remove_ltfu_10pct', 'remove_ltfu_20pct',
                               'remove_f_msm_10pct_and_msw_20pct_ltfu'
                               )

  hypothetical_interventions <<- c("high_activity_young_annual",
                                  "high_activity_young_semi_annual",
                                  "young_annual", 
                                  "young_semi_annual",
                                  "high_activity_MSM_annual",
                                  "high_activity_MSM_semi_annual",
                                  "high_activity_MSM_quarter_annual",
                                  "MSM_annual", 
                                  "MSM_semi_annual",
                                  "MSM_quarter_annual", 
                                  "high_activity_annual",
                                  "high_activity_semi_annual",
                                  "high_activity_female_young_annual",
                                  "high_activity_female_young_semi_annual",
                                  "female_young_annual",
                                  "female_young_semi_annual",
                                  '1year_high_activity_semi_annual')

  
  # interventions <<- setdiff(unique(df$interv), 'base_case')
  # saveRDS(interventions, file.path(system.file("/intervention_analysis/", package='gcRegional'), "used_interventions3.rds"))
  interventions <<- readRDS(file.path(system.file("/intervention_analysis/", package='gcRegional'), "used_interventions.rds"))
}

.render_intervention_prevalence_summary_df <- function() { 

  df <- readRDS(system.file("intervention_analysis/interventions3.rds", package='gcRegional'))

  .construct_interventions_description_df()

  int_prev_df <- df %>% 
    select(
      c(
        "site",
        "sim",
        "interv",
        "year",
        "population",
        "population_size",
        "symptomatic_cases",
        "asymptomatic_cases"
        )) %>% 
    filter(population == "All") %>% 
    select(-population) %>% 
    mutate(prevalence = (asymptomatic_cases + symptomatic_cases) / population_size) %>% 
    select(-c(population_size, asymptomatic_cases, symptomatic_cases))

  saveRDS(int_prev_df, file = file.path(
    system.file(
      package = 'gcRegional',
      "intervention_analysis/data/"
    ),
    'int_prev_df3.rds'
  ))

  int_prev_ci_df <- group_by(int_prev_df, site, interv, year) %>% 
    summarise(mean = mean(prevalence, na.rm=T),
              ci_high = quantile(prevalence, 0.975, na.rm = T),
              ci_low = quantile(prevalence, 0.025, na.rm = T))

  last_int_prev_data_points <- as.data.frame(do.call(
    rbind.data.frame,
    list(
      filter(int_prev_ci_df, site == 'BA', year == 2022),
      filter(int_prev_ci_df, site == 'SF', year == 2021)
    )))

  int_prev_ci_df$site <-
    sapply(as.character(int_prev_ci_df$site), function(x)
    c(SF = 'San Francisco', BA = 'Baltimore')[[x]])

  last_int_prev_data_points$site <-
    sapply(as.character(last_int_prev_data_points$site), function(x)
    c(SF = 'San Francisco', BA = 'Baltimore')[[x]])

  int_prev_ci_df[int_prev_ci_df$site == 'San Francisco', 'year'] <- int_prev_ci_df[int_prev_ci_df$site == 'San Francisco', 'year'] - 2016
  int_prev_ci_df[int_prev_ci_df$site == 'Baltimore', 'year'] <- int_prev_ci_df[int_prev_ci_df$site == 'Baltimore', 'year'] - 2017

  int_prev_ci_df <- merge(int_prev_ci_df, intervention_descriptions, by.x = 'interv', by.y = 'shortnames', all.x = T)

  last_int_prev_data_points[last_int_prev_data_points$site == 'San Francisco', 'year'] <-
    last_int_prev_data_points[last_int_prev_data_points$site == 'San Francisco', 'year'] - 2016
  last_int_prev_data_points[last_int_prev_data_points$site == 'Baltimore', 'year'] <-
    last_int_prev_data_points[last_int_prev_data_points$site == 'Baltimore', 'year'] - 2017

  last_int_prev_data_points <- merge(last_int_prev_data_points, intervention_descriptions, by.x = 'interv', by.y = 'shortnames', all.x = T)
  last_int_prev_data_points$names <- ordered(last_int_prev_data_points$names, levels = rev(intervention_descriptions$names))

  write.csv(x = last_int_prev_data_points, file = file.path(
    system.file(
      package = 'gcRegional',
      "intervention_analysis/data/"
    ),
    'prevalence_in_interventions_after_5_years3.csv'
  ))

  write.csv(x = int_prev_ci_df, file = file.path(
    system.file(
      package = 'gcRegional',
      "intervention_analysis/data/"
    ),
    'prevalence_in_interventions_in_all_years3.csv'
  ))

  saveRDS(last_int_prev_data_points, file = file.path(
    system.file(
      package = 'gcRegional',
      "intervention_analysis/data/"
    ),
    'last_int_prev_data_points3.rds'
  ))

  saveRDS(int_prev_ci_df, file = file.path(
    system.file(
      package = 'gcRegional',
      "intervention_analysis/data/"
    ),
    'int_prev_ci_df3.rds'
  ))
}


.render_relative_reduction_in_intervention_prevalence <- function() { 
  
  int_prev_df <- readRDS(file.path(
    system.file(
      "intervention_analysis/data/",
      package = 'gcRegional'
    ),
    'int_prev_df3.rds'
  ))

  int_prev_dcast <- dcast(int_prev_df, site + sim + year ~ interv)

  for (interv in interventions) {
      int_prev_dcast[[interv]] = int_prev_dcast[[interv]] / int_prev_dcast[['base_case']]
  }

  int_prev_dcast %<>% select(-base_case)

  relative_summary_df <- melt(int_prev_dcast, id.vars = c('site', 'sim', 'year'), variable.name = 'interv', value.name = 'prevalence')

  relative_ci_df <- relative_summary_df %>%  
    group_by(site, interv, year) %>% 
    summarise(mean = mean(prevalence, na.rm = T), 
              ci_high = quantile(prevalence, 0.975, na.rm = T),
              iq_high = quantile(prevalence, 0.75, na.rm = T),
              ci_low = quantile(prevalence, 0.025, na.rm = T),
              iq_low = quantile(prevalence, 0.25, na.rm = T))

  last_data_points_relative <- as.data.frame(do.call(
    rbind.data.frame,
    list(
      filter(relative_ci_df, site == 'BA', year == 2022),
      filter(relative_ci_df, site == 'SF', year == 2021)
    )))

  last_data_points_relative$site <-
    sapply(as.character(last_data_points_relative$site), function(x)
    c(SF = 'San Francisco', BA = 'Baltimore')[[x]])

  last_data_points_relative <- merge(last_data_points_relative, intervention_descriptions, by.x = 'interv', by.y = 'shortnames', all.x = T)
  last_data_points_relative$names <- ordered(last_data_points_relative$names, levels = rev(intervention_descriptions$names))

  write.csv(x = last_data_points_relative, file = file.path(
    system.file(
      package = 'gcRegional',
      "intervention_analysis/data/"
    ),
    'relative_prevalence_reductions_in_interventions_after_5_years3.csv'
  ))

  saveRDS(last_data_points_relative, file.path(
    system.file(
      package = 'gcRegional',
      "intervention_analysis/data/"
    ),
    'relative_prevalence_reductions_in_interventions_after_5_years3.rds'
  ))


}


.render_additional_tests_cases_averted_nnt_data <- function() { 

  .construct_interventions_description_df()

  df <- readRDS(file.path(
    system.file("intervention_analysis/", package='gcRegional'),
    "interventions3.rds"))

  # Filtering for relevant data
  additional_tests_df <- df %>% 
    dplyr::select(
      c(
        "site",
        "sim",
        "interv",
        "year",
        "population",
        "number_screened"
        )) %>% 
    filter(population == "All",
           (site == 'BA' & year == 2022) | 
             (site == 'SF' & year == 2021)
      ) %>% 
    select(-c(population, year))

  # Recast the data so that the numbers of screened individuals appear in columns specific
  # to each intervention instead of all in one column with another 'interv' column.
  additional_tests_df %<>% dcast(site + sim ~ interv, value.var = 'number_screened')

  # Compute number of additional tests relative to the base case
  for (interv in interventions) {
      additional_tests_df[[interv]] = additional_tests_df[[interv]] / additional_tests_df[['base_case']]
  }

  # Melt the data so that interv becomes a column again and 
  # so that the data is in a format useable with ggplot2.
  additional_tests_df %<>% melt(id.vars = c('site', 'sim'), variable.name = 'interv')

  # Remove 1 so that we are visualizing the relative increase (rather than relative rate).
  additional_tests_df$value <- additional_tests_df$value - 1

  saveRDS(additional_tests_df, file.path(
    system.file("intervention_analysis", package='gcRegional'),
    'additional_tests_df3.rds'))

  # Summarize the data with interquartile ranges, confidence intervals, and mean.
  additional_tests_df %<>%  
    group_by(site, interv) %>% 
    summarise(mean = mean(value, na.rm=T),
              median = median(value, na.rm = T), 
              ci_high = quantile(value, 0.975, na.rm = T),
              iq_high = quantile(value, 0.75, na.rm = T),
              ci_low = quantile(value, 0.025, na.rm = T),
              iq_low = quantile(value, 0.25, na.rm = T))


  # Select relevant cases averted data.
  cases_averted_df <- df %>% 
    select(
      c(
        "site",
        "sim",
        "interv",
        "year",
        "population",
        "incident_cases"
        )) %>% 
    filter(population == "All",
           (site == 'BA' & year == 2022) | 
             (site == 'SF' & year == 2021)
      ) %>% 
    select(-c(population, year))

  # Recast data so that each intervention's incident cases is each a column.
  cases_averted_df %<>% dcast(site + sim ~ interv, value.var = 'incident_cases')

  # Compute incident cases in interventions as relative to the base case.
  for (interv in interventions) {
      cases_averted_df[[interv]] = (cases_averted_df[['base_case']] - cases_averted_df[[interv]]) / cases_averted_df[['base_case']]
  }

  # Melt the data so that interventions' data all goes into one column again 
  # so we can use ggplot2.
  cases_averted_df %<>% melt(id.vars = c('site', 'sim'), variable.name = 'interv')

  saveRDS(cases_averted_df, file.path(
    system.file("intervention_analysis", package='gcRegional'),
    'cases_averted_df3.rds'))

  # Summarize the data with confidence intervals, interquartile ranges, and mean.
  cases_averted_df %<>%  
    group_by(site, interv) %>% 
    summarise(mean = mean(value, na.rm=T),
              median = median(value, na.rm = T), 
              ci_high = quantile(value, 0.975, na.rm = T),
              iq_high = quantile(value, 0.75, na.rm = T),
              ci_low = quantile(value, 0.025, na.rm = T),
              iq_low = quantile(value, 0.25, na.rm = T))


  #  Data needed to compute number needed to test
  nnt_df <- df %>% 
    select(
      c(
        "site",
        "sim",
        "interv",
        "year",
        "population",
        "incident_cases",
        "number_screened"
        )) %>% 
    filter(population == "All",
           (site == 'BA' & year == 2022) | 
             (site == 'SF' & year == 2021)
      ) %>% 
    select(-c(year, population))

  nnt_df_basecase <- filter(nnt_df, interv == 'base_case') %>% select(-interv)

  nnt_df <- merge(nnt_df, nnt_df_basecase, by	= c('site', 'sim'))

  # nnt_df_orig <- nnt_df 

  nnt_df %<>% mutate(
    cases_averted = incident_cases.y - incident_cases.x,
    additional_tests = number_screened.x - number_screened.y) %>% 
    select(-c(incident_cases.y, incident_cases.x, number_screened.y, number_screened.x)) %>% 
    mutate(nnt = additional_tests / cases_averted) %>% 
    filter(interv != 'base_case') %>% 
    rename(value = nnt) # %>% 
    # select(-c('cases_averted', 'additional_tests'))

  # nnt_df_2 <- nnt_df %>% mutate(
  #   cases_averted = incident_cases.y - incident_cases.x / incident_cases.y,
  #   additional_tests = (number_screened.x - number_screened.y)/number_screened.y) %>% 
  #   select(-c(incident_cases.y, incident_cases.x, number_screened.y, number_screened.x)) %>% 
  #   mutate(nnt = additional_tests / cases_averted) %>% 
  #   filter(interv != 'base_case') %>% 
  #   rename(value = nnt) # %>% 

  saveRDS(nnt_df, file.path(
    system.file("intervention_analysis", package='gcRegional'),
    'nnt_df3.rds'))

  # saveRDS(nnt_df2, file.path(
  #   system.file("intervention_analysis", package='gcRegional'),
  #   'nnt_df2.rds'))


  nnt_df %<>%  
    group_by(site, interv) %>% 
    summarise(# mean1 = mean(value, na.rm=T),
              mean = mean(additional_tests, na.rm=T) / mean(cases_averted, na.rm=T),
              median = median(value, na.rm = T), 
              ci_high = quantile(value, 0.975, na.rm = T),
              iq_high = quantile(value, 0.75, na.rm = T),
              ci_low = quantile(value, 0.025, na.rm = T),
              iq_low = quantile(value, 0.25, na.rm = T))

  # Merge the data together.
  nt_ca_df <- do.call(rbind, list(
    cbind.data.frame(variable = 'additional_tests', additional_tests_df),
    cbind.data.frame(variable = 'cases_averted', cases_averted_df),
    cbind.data.frame(variable = 'nnt', nnt_df)
  ))

  # Replace variables with human readable names.
  nt_ca_df$variable <-
    sapply(nt_ca_df$variable, function(x)
    c(additional_tests = 'Additional Tests (%)', 
      cases_averted = 'Infections Averted (%)', 
      nnt = 'Number Needed to Test')[[x]])


  # Merge intervention descriptions into our data.
  nt_ca_df <- merge(nt_ca_df, intervention_descriptions, by.x = 'interv', by.y = 'shortnames', all.x = T)

  # Order the data frame by the intervention description names.
  nt_ca_df$names <- ordered(nt_ca_df$names, rev(intervention_descriptions$names))

  # Replace 'SF' and 'BA' with their real names.
  nt_ca_df$site <-
    sapply(as.character(nt_ca_df$site), function(x)
    c(SF = 'San Francisco', BA = 'Baltimore')[[x]])

  nt_ca_df[grepl('ltfu', nt_ca_df$interv) & nt_ca_df$variable == 'Additional Tests (%)', 
           c('mean', 'ci_high', 'ci_low', 'iq_high', 'iq_low')] <-
             NA
  nt_ca_df[grepl('ltfu', nt_ca_df$interv) & nt_ca_df$variable == 'Number Needed to Test', 
           c('mean', 'ci_high', 'ci_low', 'iq_high', 'iq_low')] <-
             NA
             
  write.csv(x = nt_ca_df, file = file.path(
    system.file(
      package = 'gcRegional',
      "intervention_analysis/data/"
    ),
    'additional_tests_(pct)_cases_averted_(pct)_and_nnt3.csv'
  ))

  saveRDS(nt_ca_df, file = file.path(
    system.file(
      package = 'gcRegional',
      "intervention_analysis/data/"
    ),
    'additional_tests_(pct)_cases_averted_(pct)_and_nnt3.rds'
  ))

}
