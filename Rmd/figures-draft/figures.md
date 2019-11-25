Figures
=======

This document serves the purpose of showing where the figures available
in the manuscript are available in this package and git repository.
Additionally this document outlines the way in which each of these
figures are produced using the `gcRegional` R package.

The figures will be listed here as follows:

-   fig1 Fit to Diagnosis Rates
-   fig2 Intervention Prevalence Time Trends
-   fig3 Additional Tests and Infections Averted
-   fig4 Infections Averted in LTFU scenarios
-   calibration figures, BA and SF
-   supplement fig1 Breakdown of Modeled Incidence and Reported
    Diagnoses
-   supplement fig2 MSM Intervention Sensitivity Analysis

<!-- -->

    devtools::load_all('.')

    ## Loading gcRegional

### fig1 Fit to Diagnosis Rates

### fig2 Intervention Prevalence Time Trends

    plots <- plot_intervention_prevalence_time_trends()
    names(plots)

    ## [1] "current_strategies_interventions_prevalence"         "hypothetical_strategies_interventions_prevalence"   
    ## [3] "main_strategies_interventions_prevalence"            "hypothetical_all_activity_interventions_prevalence" 
    ## [5] "hypothetical_high_activity_interventions_prevalence"

    plots[['main_strategies_interventions_prevalence']]

![](figures_files/figure-markdown_strict/unnamed-chunk-2-1.png)

### fig3 Additional Tests and Infections Averted

##### Part A) Baltimore

    plots <- plot_number_needed_to_treat()
    plots[['additional_tests_and_cases_averted_ba_hypothetical_all_activity']]

    ## Warning: Removed 1 rows containing missing values (geom_hline).

![](figures_files/figure-markdown_strict/unnamed-chunk-3-1.png)

##### Part B) San Francisco

    plots[['additional_tests_and_cases_averted_sf_hypothetical_all_activity']]

    ## Warning: Removed 1 rows containing missing values (geom_hline).

![](figures_files/figure-markdown_strict/unnamed-chunk-4-1.png)

### fig4 Infections Averted in LTFU scenarios

    plots[['cases_averted_ltfu']]

![](figures_files/figure-markdown_strict/unnamed-chunk-5-1.png)

### supplement fig1 Breakdown of Modeled Incidence and Reported Diagnoses

    plot_comparison_of_model_incidence_and_diagnoses()

    ## Warning: Removed 14 rows containing missing values (geom_errorbar).

![](figures_files/figure-markdown_strict/unnamed-chunk-6-1.png)

### supplement fig2 MSM Intervention Sensitivity Analysis

    plots <- plot_msm_interventions()
    plots[['ba']]

    ## Warning: Removed 50 rows containing missing values (geom_path).

![](figures_files/figure-markdown_strict/unnamed-chunk-7-1.png)

    plots[['sf']]

    ## Warning: Removed 50 rows containing missing values (geom_path).

![](figures_files/figure-markdown_strict/unnamed-chunk-7-2.png)

    interventions_msm <- readRDS(file.path(
      system.file(
        'intervention_analysis/',
        package = 'gcRegional'), 'interventions_msm.rds'))

    ggplot(filter(interventions_msm, site == 'BA', year == 2023), aes(x=additional_tests, y=cases_averted)) + 
      geom_point(alpha=0.5) + 
      theme_minimal() + 
      ggtitle("Baltimore")

    ## Warning: Removed 1000 rows containing missing values (geom_point).

![](figures_files/figure-markdown_strict/unnamed-chunk-8-1.png)

    ggplot(filter(interventions_msm, site == 'SF', year == 2022), aes(x=additional_tests, y=cases_averted)) + 
      geom_point(alpha=0.5) + 
      theme_minimal() + 
      ggtitle("San Francisco")

    ## Warning: Removed 1068 rows containing missing values (geom_point).

![](figures_files/figure-markdown_strict/unnamed-chunk-8-2.png)

    ggplot(filter(interventions_msm, site == 'BA', year == 2023), aes(x=additional_tests, y=cases_averted)) + 
      geom_point(alpha=0.5) + 
      theme_minimal() + 
      geom_smooth() + 
      ggtitle("Baltimore")

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 1000 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 1000 rows containing missing values (geom_point).

![](figures_files/figure-markdown_strict/unnamed-chunk-8-3.png)

    ggplot(filter(interventions_msm, site == 'SF', year == 2022), aes(x=additional_tests, y=cases_averted)) + 
      geom_point(alpha=0.5) + 
      theme_minimal() + 
      geom_smooth() + 
      ggtitle("San Francisco")

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 1068 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 1068 rows containing missing values (geom_point).

![](figures_files/figure-markdown_strict/unnamed-chunk-8-4.png)
