
#' Plot Comparison of Modeled Incidence vs. Reported and Diagnosed Infections
#'
#' The figure this code generates corresponds to supplement 3 figure 1. 
#' This code relies on the intervention analysis data having been properly 
#' formatted and stored by the function `.render_incidence_vs_reported_cases_data()`
#' 
#' @param save Indicates whether or not the plot should be ggsaved into the package.
#' 
#' 
plot_comparison_of_model_incidence_and_diagnoses <- function(save = F) { 
  diag_inc_comparison_df <- readRDS(file.path(
    system.file("intervention_analysis/data/", package='gcRegional'),
    "diag_inc_comparison_df.rds"))

  diag_inc_comparison_df_errorbars <- readRDS(file.path(
    system.file("intervention_analysis/data/", package='gcRegional'),
    "diag_inc_comparison_df_errorbars.rds"))

  model_estimated_incident_infections_breakdown <- ggplot(
   diag_inc_comparison_df,
    aes(
    x = population,
    y = 100*mean,
    ymax = 100*ci_high,
    ymin = 100*ci_low,
    lower = 100*iq_low,
    upper = 100*iq_high,
    fill = variable,
    group = variable
    )) + 
    geom_bar(stat = 'identity', position=position_dodge()) + 
    geom_errorbar(data = diag_inc_comparison_df_errorbars,
                  stat = 'identity', position=position_dodge(.85), width = 0.2, alpha=0.5) + 
    facet_wrap(~site) + 
    ylab("Percentage") + 
    xlab("Demographic Groups") + 
    ggtitle("Model Estimated Proportional Breakdown of Incident Infections", subtitle='Compared to Diagnosed and Reported Cases') + 
    theme_bw() + 
    ylim(c(0,100)) + 
    guides(fill=guide_legend(title='')) +
    scale_fill_manual(values = c('model_estimate' = '#2b8cbe', 'reported_cases' = '#fdbb84'),
                      labels = c('model_estimate' = 'Model Estimated Proportional\nBreakdown of Incident Infections',
                           'reported_cases' = 'Proportional Breakdown of \nDiagnosed Cases')) + 
    scale_color_manual(values=c('target'='red', labels='Breakdown of Reported Cases')) + 
    scale_shape_manual(values=c('target'=15), labels=c('Breakdown of Reported Cases')) + 
    theme(axis.text.x = element_text(angle = 25, vjust = 0.5), legend.position = 'bottom') 




  if (save) { 
    ggsave(plot = model_estimated_incident_infections_breakdown,
           filename = file.path(
              system.file("intervention_analysis/figures/", package = 'gcRegional'),
              'distribution_of_cases.png'
              ),
           height = unit(5, 'in'))
  }

  return(model_estimated_incident_infections_breakdown)
}


#' Plot Time Trends of Prevalence in Interventions during the Intervention Time Period
#' 
#' 
plot_intervention_prevalence_time_trends <- function(save = F, bw = F) { 

  require(ggrepel)
  
  .construct_interventions_description_df()

  last_int_prev_data_points <- 
  readRDS(file.path(
    system.file(
      package = 'gcRegional',
      "intervention_analysis/data/"
    ),
    'last_int_prev_data_points3.rds'
  ))

  int_prev_ci_df <- readRDS(file.path(
    system.file(
      package = 'gcRegional',
      "intervention_analysis/data/"
    ),
    'int_prev_ci_df3.rds'
  ))

  make_interventions_prevalence_plot <- function(interventions_list) {
    ggplot(
      filter(
      int_prev_ci_df,
      year >= 0,
      year <= 5,
      interv %in% interventions_list
      ),
      aes(x = year, y = 100 * mean, color = names)
      ) +
      geom_line() +
    facet_wrap(~site) + 
    scale_x_continuous(
      breaks = seq(0, 5, by = 1), 
      limits = c(0, 9)) +
    geom_label_repel(
      data = filter(last_int_prev_data_points, interv %in% interventions_list),
      aes(
      x = year,
      y = 100*mean,
      group = interv,
      label = names,
      # color = names,
      fill = names
      ),
      color = 'black',
      segment.color = 'black',
      segment.size = 0.1,
      segment.alpha = 0.8,
      # fontface = 'bold',
      alpha = 0.95,
      nudge_x = 2.5,
      direction = 'y',
      size = 2.5
      ) +
    expand_limits(y = 0) + 
    ylab("Prevalence (%)") +
    xlab("Years After Intervention Start") + 
    ggtitle("Prevalence Estimates During the Intervention Time Period") + 
    theme_bw() + 
    theme(legend.position = 'none', panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank()) +
    # if (bw) { 
    #   scale_color_manual(values = setNames(rep('black', nrow(intervention_descriptions)), intervention_descriptions$names))
    # } else { 
      scale_color_manual(values = setNames(intervention_descriptions$color, intervention_descriptions$names)) + 
    # } +  
    # if (bw) { 
    #   scale_fill_manual(values = '#FFFFFF')
    # } else {
      scale_fill_manual(values = setNames(sapply(intervention_descriptions$color,
                                                 colorspace::lighten, amount=.85),
                                          intervention_descriptions$names))
    # }
  }

  plots <- list()

  plots[['current_strategies_interventions_prevalence']] <- make_interventions_prevalence_plot(c(empirical_interventions, 'base_case'))

  plots[['hypothetical_strategies_interventions_prevalence']] <- make_interventions_prevalence_plot(c('base_case', hypothetical_interventions))

  main_strategies <- c('base_case', empirical_interventions,
                       hypothetical_interventions[!grepl('high', hypothetical_interventions)],
                       '1year_high_activity_semi_annual')
                       
  plots[['main_strategies_interventions_prevalence']] <- make_interventions_prevalence_plot(main_strategies)

  hypothetical_ints_all_activity <- hypothetical_interventions[!grepl("high", hypothetical_interventions)]
  plots[['hypothetical_all_activity_interventions_prevalence']] <- make_interventions_prevalence_plot(c('base_case', hypothetical_ints_all_activity, '1year_high_activity_semi_annual'))


  hypothetical_ints_high_activity <- hypothetical_interventions[grepl("high", hypothetical_interventions)]
  plots[['hypothetical_high_activity_interventions_prevalence']] <- make_interventions_prevalence_plot(c('base_case', setdiff(hypothetical_ints_high_activity, c('1year_high_activity_semi_annual', 'high_activity_annual', 'high_activity_semi_annual'))))

  if (save) { 

    bw_suffix <- if (bw) '_bw' else ''

    ggsave(file.path(
      system.file("intervention_analysis/figures/", package = 'gcRegional'),
      paste0('intervention_prevalence_current_strategies', bw_suffix, '.png')
      ),
      height = unit(5, 'in'),
      width = unit(9, 'in'),
      plot = plots[['current_strategies_interventions_prevalence']])

    ggsave(file.path(
      system.file("intervention_analysis/figures/", package = 'gcRegional'),
      paste0('intervention_prevalence_hypothetical', bw_suffix, '.png')
      ),
      height = unit(5, 'in'),
      width = unit(9, 'in'),
      plot = plots[['hypothetical_strategies_interventions_prevalence']])

    ggsave(file.path(
      system.file("intervention_analysis/figures/", package = 'gcRegional'),
      paste0('intervention_prevalence_main', bw_suffix, '.png')
      ),
      height = unit(7, 'in'),
      width = unit(9, 'in'),
      plot = plots[['main_strategies_interventions_prevalence']])

    ggsave(file.path(
      system.file("intervention_analysis/figures/", package = 'gcRegional'),
      paste0('intervention_prevalence_hypothetical_all_activity', bw_suffix, '.png')
      ),
      height = unit(5, 'in'),
      width = unit(9, 'in'),
      plot = plots[['hypothetical_all_activity_interventions_prevalence']])

    ggsave(file.path(
      system.file("intervention_analysis/figures/", package = 'gcRegional'),
      paste0('intervention_prevalence_hypothetical_high_activity', bw_suffix, '.png')
      ),
      height = unit(5, 'in'),
      width = unit(10, 'in'),
      plot = plots[['hypothetical_high_activity_interventions_prevalence']])
  } 

  return(plots)
}

#' Plot Relative Reductions in Prevalence in Interventions
plot_intervention_relative_reductions_in_prevalence <- function(save=F) { 
  .construct_interventions_description_df()

  last_data_points_relative <- readRDS(file.path(
    system.file(
      package = 'gcRegional',
      "intervention_analysis/data/"
    ),
    'relative_prevalence_reductions_in_interventions_after_5_years.rds'
  ))

  make_relative_prevalence_reduction_plot <- function(interventions_list) {
    ggplot(
      filter(last_data_points_relative, interv %in% interventions_list),
      aes(
      x = names,
      middle = 100 - 100 * mean,
      ymax = 100 - 100 * ci_high,
      ymin = 100 - 100 * ci_low,
      upper = 100 - 100 * iq_high,
      lower = 100 - 100 * iq_low,
      fill = names
      )
      ) +
      geom_boxplot(stat = 'identity') +
      facet_wrap( ~ site) +
      coord_flip() +
      theme_bw() +
      theme(legend.position = 'none',
      axis.title.x = element_text(hjust = 1), 
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
      ) +
      xlab("") +
      ylab("Prevalence Reduction as a Percentage of Base Case Prevalence") +
      ylim(y = c(0, 100)) +
      ggtitle("Relative Reduction in Prevalence") +
      scale_fill_manual(values = setNames(
      sapply(intervention_descriptions$color, colorspace::lighten),
      intervention_descriptions$names
      ))
      
  }

  plots <- list()

  plots[['empirical_relative_reductions']] <- make_relative_prevalence_reduction_plot(empirical_interventions)

  plots[['hypothetical_relative_reductions ']]<- make_relative_prevalence_reduction_plot(hypothetical_interventions)

  plots[['hypothetical_all_activity_relative_reductions ']]<- make_relative_prevalence_reduction_plot(
    c(hypothetical_interventions[!grepl('high', hypothetical_interventions)], '1year_high_activity_semi_annual'))

  plots[['hypothetical_high_activity_relative_reductions ']]<-
    make_relative_prevalence_reduction_plot(setdiff(hypothetical_interventions[grepl('high',
      hypothetical_interventions)], c('1year_high_activity_semi_annual',
      'high_activity_annual', 'high_activity_semi_annual')))

  if (save) { 
    ggsave(file.path(
      system.file("intervention_analysis/figures/", package = 'gcRegional'),
      'prevalence_ratio_in_last_year_empirical.png'
      ),
      height = unit(5, 'in'), 
      width = unit(9, 'in'),
      plot = plots[['empirical_relative_reductions']])

    ggsave(file.path(
      system.file("intervention_analysis/figures/", package = 'gcRegional'),
      'prevalence_ratio_in_last_year_hypothetical.png'
      ),
      height = unit(5, 'in'), 
      width = unit(9, 'in'),
      plot = plots[['hypothetical_relative_reductions']])

  ggsave(file.path(
    system.file("intervention_analysis/figures/", package = 'gcRegional'),
    'prevalence_ratio_in_last_year_hypothetical_all_activity.png'
    ),
    height = unit(5, 'in'), 
    width = unit(9, 'in'),
    plot = plots[['hypothetical_all_activity_relative_reductions']])

  ggsave(file.path(
    system.file("intervention_analysis/figures/", package = 'gcRegional'),
    'prevalence_ratio_in_last_year_hypothetical_high_activity.png'
    ),
    height = unit(5, 'in'), 
    width = unit(9, 'in'),
    plot = plots[['hypothetical_high_activity_relative_reductions']])
  }

  return(plots)
}


plot_number_needed_to_treat <- function(save=F) { 

  .construct_interventions_description_df()

  nt_ca_df <- readRDS(file.path(
    system.file(
      package = 'gcRegional',
      "intervention_analysis/data/"
    ),
    'additional_tests_(pct)_cases_averted_(pct)_and_nnt3.rds'
  ))

  make_additional_tests_and_cases_averted_plot <- function(site, interventions_list, subtitle, include_nnt = F) {
    ggplot(filter(nt_ca_df, 
                  site == !! site, 
                  variable %in% c('Additional Tests (%)', 
                                  'Infections Averted (%)', 
                                  if (include_nnt) 'Number Needed to Test' else NULL), 
                  interv %in% interventions_list)) + 
    geom_hline(data = data.frame(variable = c('Additional Tests (%)', 'Infections Averted (%)'), yintercept = c(NA, 0)), 
               mapping = aes(yintercept=yintercept), color='red', alpha=0.5) + 
    geom_boxplot(
      stat = 'identity',
      aes(
      x = names,
      ymax = 100*ci_high,
      upper = 100*iq_high,
      middle = 100*median,
      lower = 100*iq_low,
      ymin = 100*ci_low,
      fill = names
      )
      ) + 
    ggtitle(paste0(c(`Baltimore` = 'A) ', `San Francisco` = 'B) ')[[site]], "Intervention Outcomes in ", site)) + 
    theme_bw() + 
    facet_wrap(~variable, scales = 'free_x') + 
    coord_flip()  +
    xlab("") + 
    expand_limits(y=0) + 
    theme(axis.title.x = element_text(hjust = 1), legend.position = 'none', 
          axis.text.x = element_text(hjust=.5),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          plot.margin = margin(10, 15, 10, 10)) + 
    ylab("Additional Tests and Infections Averted as a Percentage of the Base Case") + 
      scale_fill_manual(values = setNames(
        sapply(intervention_descriptions$color, colorspace::lighten),
        intervention_descriptions$names
      ))
  }

  plots <- list()

  # Additional Tests and Cases Averted - Baltimore, Empirical Interventions
  plots[['additional_tests_and_cases_averted_ba_empirical']] <-
    make_additional_tests_and_cases_averted_plot('Baltimore', empirical_interventions)

  # Cases Averted in the LTFU scenarios
  plots[['cases_averted_ltfu']] <-  
  ggplot(filter(nt_ca_df, 
                variable == 'Infections Averted (%)', 
                interv %in% empirical_interventions)) + 
    geom_boxplot(
      stat = 'identity',
      aes(
      x = names,
      ymax = 100*ci_high,
      upper = 100*iq_high,
      middle = 100*median,
      lower = 100*iq_low,
      ymin = 100*ci_low,
      fill = names
      )
      ) + 
    ggtitle("Infections Averted in LTFU Interventions") + 
    theme_bw() + 
    facet_wrap(~site) + 
    coord_flip()  +
    xlab("") + 
    expand_limits(y=0) + 
    theme(axis.title.x = element_text(hjust = 1), legend.position = 'none', 
          axis.text.x = element_text(hjust=.5),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          plot.margin = margin(10, 15, 10, 10)) + 
    ylab("Infections Averted as a Percentage of the Base Case") + 
      scale_fill_manual(values = setNames(
        sapply(intervention_descriptions$color, colorspace::lighten),
        intervention_descriptions$names
      ))

  # Additional Tests and Cases Averted - San Francisco, Empirical Interventions
  additional_tests_and_cases_averted_sf_empirical <-
    make_additional_tests_and_cases_averted_plot('San Francisco', empirical_interventions)

  # Additional Tests and Cases Averted - Baltimore, Hypothetical Interventions
  plots[['additional_tests_and_cases_averted_ba_hypothetical']] <-
    make_additional_tests_and_cases_averted_plot('Baltimore', hypothetical_interventions)

  # Additional Tests and Cases Averted - San Francisco, Hypothetical Interventions
  plots[['additional_tests_and_cases_averted_sf_hypothetical']] <-
    make_additional_tests_and_cases_averted_plot('San Francisco', hypothetical_interventions)

  # Additional Tests and Cases Averted - Baltimore, Hypothetical Interventions without High Activity Targeted
  plots[['additional_tests_and_cases_averted_ba_hypothetical_all_activity']] <-
    make_additional_tests_and_cases_averted_plot('Baltimore',
    c(hypothetical_interventions[!grepl('high', hypothetical_interventions)],
    '1year_high_activity_semi_annual'), subtitle = 'Hypothetical
    Interventions')

  # Additional Tests and Cases Averted - San Francisco, Hypothetical Interventions without High Activity Targeted
  plots[['additional_tests_and_cases_averted_sf_hypothetical_all_activity']] <-
    make_additional_tests_and_cases_averted_plot('San Francisco',
    c(hypothetical_interventions[!grepl('high', hypothetical_interventions)],
    '1year_high_activity_semi_annual'), subtitle = 'Hypothetical
    Interventions')

  # Additional Tests and Cases Averted - Baltimore, Hypothetical Interventions - Only High Activity Targeted
  plots[['additional_tests_and_cases_averted_ba_hypothetical_high_activity']] <-
    make_additional_tests_and_cases_averted_plot('Baltimore',
    setdiff(hypothetical_interventions[grepl('high',
    hypothetical_interventions)], c('1year_high_activity_semi_annual',
    'high_activity_annual', 'high_activity_semi_annual')), subtitle = 'Risk
    Targeted Hypothetical Interventions')

  # Additional Tests and Cases Averted - San Francisco, Hypothetical Interventions - Only High Activity Targeted
  plots[['additional_tests_and_cases_averted_sf_hypothetical_high_activity']] <-
    make_additional_tests_and_cases_averted_plot('San Francisco',
    setdiff(hypothetical_interventions[grepl('high',
    hypothetical_interventions)], c('1year_high_activity_semi_annual',
    'high_activity_annual', 'high_activity_semi_annual')), subtitle = 'Risk
    Targeted Hypothetical Interventions')


  if (save) {
    ggsave(
      file.path(
      system.file("intervention_analysis/figures/", package = 'gcRegional'),
      'number_tested_cases_averted_nnt_baltimore_empirical.png'
      ),
      height = unit(6, 'in'),
      width = unit(9, 'in'),
      plot = plots[['additional_tests_and_cases_averted_ba_empirical']] )

    ggsave(
      file.path(
      system.file("intervention_analysis/figures/", package = 'gcRegional'),
      'cases_averted_ltfu_scenarios.png'
      ),
      height = unit(6, 'in'),
      width = unit(9, 'in'),
      plot = plots[['cases_averted_ltfu']] )
      
    ggsave(
      file.path(
      system.file("intervention_analysis/figures/", package = 'gcRegional'),
      'number_tested_and_cases_averted_sanfrancisco_empirical.png'
      ),
      height = unit(6, 'in'),
      width = unit(9, 'in'),
      plot = additional_tests_and_cases_averted_sf_empirical
      )

    ggsave(
      file.path(
      system.file("intervention_analysis/figures/", package = 'gcRegional'),
      'number_tested_and_cases_averted_baltimore_hypothetical.png'
      ),
      height = unit(6, 'in'),
      width = unit(9, 'in'),
      plot = plots[['additional_tests_and_cases_averted_ba_hypothetical
     ']] )

    ggsave(
      file.path(
      system.file("intervention_analysis/figures/", package = 'gcRegional'),
      'number_tested_and_cases_averted_sanfrancisco_hypothetical.png'
      ),
      height = unit(6, 'in'),
      width = unit(9, 'in'),
      plot = plots[['additional_tests_and_cases_averted_sf_hypothetical
     ']] )

    ggsave(
      file.path(
      system.file("intervention_analysis/figures/", package = 'gcRegional'),
      'number_tested_and_cases_averted_baltimore_hypothetical_all_activity.png'
      ),
      height = unit(6, 'in'),
      width = unit(9, 'in'),
      plot = plots[['additional_tests_and_cases_averted_ba_hypothetical_all_activity']] )

    ggsave(
      file.path(
      system.file("intervention_analysis/figures/", package = 'gcRegional'),
      'number_tested_and_cases_averted_sanfrancisco_hypothetical_all_activity.png'
      ),
      height = unit(6, 'in'),
      width = unit(9, 'in'),
      plot = plots[['additional_tests_and_cases_averted_sf_hypothetical_all_activity']] )


    ggsave(
      file.path(
      system.file("intervention_analysis/figures/", package = 'gcRegional'),
      'number_tested_and_cases_averted_baltimore_hypothetical_high_activity.png'
      ),
      height = unit(6, 'in'),
      width = unit(9, 'in'),
      plot = plots[['additional_tests_and_cases_averted_ba_hypothetical_high_activity']] )

    ggsave(
      file.path(
      system.file("intervention_analysis/figures/", package = 'gcRegional'),
      'number_tested_and_cases_averted_sanfrancisco_hypothetical_high_activity.png'
      ),
      height = unit(6, 'in'),
      width = unit(9, 'in'),
      plot = plots[['additional_tests_and_cases_averted_sf_hypothetical_high_activity']] )

  }
  return(plots)
}


#' 
#' Second Version with Jitter of Outcome Boxplots

plot_number_needed_to_treat2 <- function(save=F) { 

  .construct_interventions_description_df()

  cases_averted_df <- readRDS(file.path(
    system.file("intervention_analysis", package='gcRegional'),
    'cases_averted_df.rds'))

  additional_tests_df <- readRDS(file.path(
    system.file("intervention_analysis", package='gcRegional'),
    'additional_tests_df.rds'))

  # nnt_df <- readRDS(file.path(
  #   system.file('intervention_analysis/nnt_df2.rds', package='gcRegional')))

  nnt_df <- rbind.data.frame(
    cbind.data.frame(cases_averted_df, variable = 'Infections Averted (%)'),
    cbind.data.frame(additional_tests_df, variable = 'Additional Tests (%)'))

  nnt_df %<>% merge(intervention_descriptions, by.x = "interv", by.y = "shortnames")

  nnt_df$variable <- factor(nnt_df$variable, levels = c('Infections Averted (%)', 'Additional Tests (%)'))

  nnt_df$interv <- factor(nnt_df$interv, levels = unique(intervention_descriptions$shortnames))
  nnt_df$names <- factor(nnt_df$names, levels = rev(unique(intervention_descriptions$names)))

  # Helper Function for Plot Making
  make_additional_tests_and_cases_averted_plot <- function(
    city = c("Baltimore", "San Francisco")[rbinom(1,1,.5)+1], 
    interventions_list = hypothetical_interventions, 
    subtitle = NULL, 
    include_nnt = F) {

    require(rlang)
    city <- rlang::enquo(city)

    df <- dplyr::filter(nnt_df, 
                  site == !! c(`Baltimore` = 'BA', `San Francisco`='SF')[[quo_name(city)]], 
                  interv %in% interventions_list)

    # df <- reshape2::melt(df, id.vars = c('site', 'sim', 'interv'))
    
    ggplot(df, aes(x=names, y = value*100)) + 
    geom_hline(data = data.frame(variable = c('Additional Tests (%)', 'Infections Averted (%)'), yintercept = c(NA, 0)), 
               mapping = aes(yintercept=yintercept), color='red', alpha=0.5) + 
    geom_jitter(
    data = df %>% filter(
      sim %in% sample.int(n=1000, size=250)),
    mapping = aes(color = names), width = 0.25, alpha = 0.5) + 
    geom_boxplot(
      width = .5, 
      outlier.shape = NA, 
      alpha = 0.5
      ) + 
    ggtitle(paste0(c(`Baltimore` = 'A) ', `San Francisco` = 'B) ')[[rlang::quo_name(city)]], "Intervention Outcomes in ", rlang::quo_name(city))) + 
    theme_bw() + 
    facet_wrap(~variable, scales = 'free_x') + 
    coord_flip()  +
    xlab("") + 
    expand_limits(y=0) + 
    theme(axis.title.x = element_text(hjust = 1), legend.position = 'none', 
          axis.text.x = element_text(hjust=.5),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          plot.margin = margin(10, 15, 10, 10)) + 
    ylab("Additional Tests and Infections Averted as a Percentage of the Base Case") + 
    scale_color_manual(values = setNames(
      sapply(
      intervention_descriptions$color, colorspace::lighten, .5),
      intervention_descriptions$names
    ))
  }

  plots <- list()

  # Additional Tests and Cases Averted - Baltimore, Empirical Interventions
  # plots[['additional_tests_and_cases_averted_ba_empirical']] <-
  #   make_additional_tests_and_cases_averted_plot('Baltimore', empirical_interventions)

  # Cases Averted in the LTFU scenarios
  plots[['cases_averted_ltfu']] <-  
    ggplot(
      data = filter(nnt_df, 
                  variable == 'Infections Averted (%)', 
                  interv %in% empirical_interventions) %>% 
      mutate(site = recode(site, BA = 'Baltimore', SF = 'San Francisco')),
      mapping = aes(x=names, y = value*100)) + 
    geom_hline(data = data.frame(variable = c('Infections Averted (%)'), yintercept = c(NA, 0)), 
               mapping = aes(yintercept=yintercept), color='red', alpha=0.5) + 
    geom_jitter(
    data = nnt_df %>% filter(
      sim %in% sample.int(n=1000, size=250),
      variable == 'Infections Averted (%)', 
      interv %in% empirical_interventions) %>% 
      mutate(site = recode(site, BA = 'Baltimore', SF = 'San Francisco')),
    mapping = aes(color = names), width = 0.25, alpha = 0.5) + 
    geom_boxplot(
      width = .5, 
      outlier.shape = NA, 
      alpha = 0.5
      ) + 
    ggtitle("Infections Averted in LTFU Interventions") + 
    theme_bw() + 
    facet_wrap(~site) + 
    coord_flip()  +
    xlab("") + 
    expand_limits(y=0) + 
    theme(axis.title.x = element_text(hjust = 1), legend.position = 'none', 
          axis.text.x = element_text(hjust=.5),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          plot.margin = margin(10, 15, 10, 10)) + 
    ylab("Infections Averted as a Percentage of the Base Case") + 
    scale_color_manual(values = setNames(
      sapply(intervention_descriptions$color, colorspace::lighten, .5),
      intervention_descriptions$names
    ))

  # ggplot() + 
  #   geom_boxplot(
  #     stat = 'identity',
  #     aes(
  #     x = names,
  #     ymax = 100*ci_high,
  #     upper = 100*iq_high,
  #     middle = 100*median,
  #     lower = 100*iq_low,
  #     ymin = 100*ci_low,
  #     fill = names
  #     )
  #     ) + 
  #   theme_bw() + 
  #   facet_wrap(~site) + 
  #   coord_flip()  +
  #   xlab("") + 
  #   expand_limits(y=0) + 
  #   theme(axis.title.x = element_text(hjust = 1), legend.position = 'none', 
  #         axis.text.x = element_text(hjust=.5),
  #         panel.grid.major.y = element_blank(),
  #         panel.grid.minor.y = element_blank(),
  #         plot.margin = margin(10, 15, 10, 10)) + 
  #   ylab("Infections Averted as a Percentage of the Base Case") + 
  #     scale_fill_manual(values = setNames(
  #       sapply(intervention_descriptions$color, colorspace::lighten),
  #       intervention_descriptions$names
  #     ))

  # Additional Tests and Cases Averted - San Francisco, Empirical Interventions
  # additional_tests_and_cases_averted_sf_empirical <-
  #   make_additional_tests_and_cases_averted_plot('San Francisco', empirical_interventions)

  # Additional Tests and Cases Averted - Baltimore, Hypothetical Interventions
  plots[['additional_tests_and_cases_averted_ba_hypothetical']] <-
    make_additional_tests_and_cases_averted_plot('Baltimore', hypothetical_interventions)

  # Additional Tests and Cases Averted - San Francisco, Hypothetical Interventions
  plots[['additional_tests_and_cases_averted_sf_hypothetical']] <-
    make_additional_tests_and_cases_averted_plot('San Francisco', hypothetical_interventions)

  # Additional Tests and Cases Averted - Baltimore, Hypothetical Interventions without High Activity Targeted
  plots[['additional_tests_and_cases_averted_ba_hypothetical_all_activity']] <-
    make_additional_tests_and_cases_averted_plot('Baltimore',
    c(hypothetical_interventions[!grepl('high', hypothetical_interventions)],
    '1year_high_activity_semi_annual'), subtitle = 'Hypothetical
    Interventions')

  # Additional Tests and Cases Averted - San Francisco, Hypothetical Interventions without High Activity Targeted
  plots[['additional_tests_and_cases_averted_sf_hypothetical_all_activity']] <-
    make_additional_tests_and_cases_averted_plot('San Francisco',
    c(hypothetical_interventions[!grepl('high', hypothetical_interventions)],
    '1year_high_activity_semi_annual'), subtitle = 'Hypothetical
    Interventions')

  # Additional Tests and Cases Averted - Baltimore, Hypothetical Interventions - Only High Activity Targeted
  plots[['additional_tests_and_cases_averted_ba_hypothetical_high_activity']] <-
    make_additional_tests_and_cases_averted_plot('Baltimore',
    setdiff(hypothetical_interventions[grepl('high',
    hypothetical_interventions)], c('1year_high_activity_semi_annual',
    'high_activity_annual', 'high_activity_semi_annual')), subtitle = 'Risk
    Targeted Hypothetical Interventions')

  # Additional Tests and Cases Averted - San Francisco, Hypothetical Interventions - Only High Activity Targeted
  plots[['additional_tests_and_cases_averted_sf_hypothetical_high_activity']] <-
    make_additional_tests_and_cases_averted_plot('San Francisco',
    setdiff(hypothetical_interventions[grepl('high',
    hypothetical_interventions)], c('1year_high_activity_semi_annual',
    'high_activity_annual', 'high_activity_semi_annual')), subtitle = 'Risk
    Targeted Hypothetical Interventions')


  if (save) {
    # ggsave(
    #   file.path(
    #   system.file("intervention_analysis/figures/", package = 'gcRegional'),
    #   'number_tested_cases_averted_nnt_baltimore_empirical.png'
    #   ),
    #   height = unit(6, 'in'),
    #   width = unit(9, 'in'),
    #   plot = plots[['additional_tests_and_cases_averted_ba_empirical']] )

    ggsave(
      file.path(
      system.file("intervention_analysis/figures/", package = 'gcRegional'),
      'cases_averted_ltfu_scenarios.png'
      ),
      height = unit(6, 'in'),
      width = unit(9, 'in'),
      plot = plots[['cases_averted_ltfu']] )
      
    # ggsave(
    #   file.path(
    #   system.file("intervention_analysis/figures/", package = 'gcRegional'),
    #   'number_tested_and_cases_averted_sanfrancisco_empirical.png'
    #   ),
    #   height = unit(6, 'in'),
    #   width = unit(9, 'in'),
    #   plot = additional_tests_and_cases_averted_sf_empirical
    #   )

    ggsave(
      file.path(
      system.file("intervention_analysis/figures/", package = 'gcRegional'),
      'number_tested_and_cases_averted_baltimore_hypothetical.png'
      ),
      height = unit(6, 'in'),
      width = unit(9, 'in'),
      plot = plots[['additional_tests_and_cases_averted_ba_hypothetical
     ']] )

    ggsave(
      file.path(
      system.file("intervention_analysis/figures/", package = 'gcRegional'),
      'number_tested_and_cases_averted_sanfrancisco_hypothetical.png'
      ),
      height = unit(6, 'in'),
      width = unit(9, 'in'),
      plot = plots[['additional_tests_and_cases_averted_sf_hypothetical
     ']] )

    ggsave(
      file.path(
      system.file("intervention_analysis/figures/", package = 'gcRegional'),
      'number_tested_and_cases_averted_baltimore_hypothetical_all_activity.png'
      ),
      height = unit(6, 'in'),
      width = unit(9, 'in'),
      plot = plots[['additional_tests_and_cases_averted_ba_hypothetical_all_activity']] )

    ggsave(
      file.path(
      system.file("intervention_analysis/figures/", package = 'gcRegional'),
      'number_tested_and_cases_averted_sanfrancisco_hypothetical_all_activity.png'
      ),
      height = unit(6, 'in'),
      width = unit(9, 'in'),
      plot = plots[['additional_tests_and_cases_averted_sf_hypothetical_all_activity']] )


    ggsave(
      file.path(
      system.file("intervention_analysis/figures/", package = 'gcRegional'),
      'number_tested_and_cases_averted_baltimore_hypothetical_high_activity.png'
      ),
      height = unit(6, 'in'),
      width = unit(9, 'in'),
      plot = plots[['additional_tests_and_cases_averted_ba_hypothetical_high_activity']] )

    ggsave(
      file.path(
      system.file("intervention_analysis/figures/", package = 'gcRegional'),
      'number_tested_and_cases_averted_sanfrancisco_hypothetical_high_activity.png'
      ),
      height = unit(6, 'in'),
      width = unit(9, 'in'),
      plot = plots[['additional_tests_and_cases_averted_sf_hypothetical_high_activity']] )

  }
  return(plots)
}
