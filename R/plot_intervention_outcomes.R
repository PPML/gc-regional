#' Plot the Intervention Simulation Outcomes
plot_interventions <- function() {
  require(ggplot2)
  require(here)
  require(dplyr)
  require(ggrepel)

  df <- readRDS(here("output/09-10-18/interventions_df.rds"))

  # Prevalence
  prev_df <- select(df, site, interv, years_out, prevalence) %>%
    group_by(site, years_out, interv) %>%
    summarize(
      mean = mean(prevalence, na.rm = T),
      ci_high = quantile(prevalence, .975, na.rm = T),
      ci_low = quantile(prevalence, .025, na.rm = T)
    )

  sf_prevalence <-
    ggplot(
      filter(prev_df, site == 'SF'),
      aes(x = years_out + 2016,
          y = mean,
          ymax = ci_high,
          ymin = ci_low,
          fill = interv,
          color = interv)) +
      # geom_ribbon(alpha = 0.1, size = 0) +
      geom_line() +
      geom_text_repel(
        data =
          filter(prev_df, site == 'SF', years_out == 5),
        mapping =
          aes(x = 2016 + years_out + .05, y = mean, label = interv),
        size = 3,
        nudge_x = .7,
        direction = "y",
        segment.size = .2,
        hjust = 0
      ) +
      ggtitle("Prevalence in San Francisco", subtitle = "Simulating 5 years of hypothetical interventions") +
      ylab("Prevalence (0-1 scale)") +
      xlab("Year") +
      coord_cartesian(xlim = c(2017, 2025), ylim = c(0.006, 0.012)) +
      theme(legend.position = 'none')

  ba_prevalence <-
    ggplot(
      filter(prev_df, site == 'BA'),
      aes(x = years_out + 2017,
          y = mean,
          ymax = ci_high,
          ymin = ci_low,
          fill = interv,
          color = interv)) +
    # geom_ribbon(alpha = 0.1, size = 0) +
    geom_line() +
    geom_text_repel(
      data =
        filter(prev_df, site == 'BA', years_out == 5),
      mapping =
        aes(x = 2017 + years_out + .05, y = mean, label = interv),
      size = 3,
      nudge_x = .7,
      direction = "y",
      segment.size = .2,
      hjust = 0
    ) +
    ggtitle("Prevalence in Baltimore", subtitle = "Simulating 5 years of hypothetical interventions") +
    ylab("Prevalence (0-1 scale)") +
    xlab("Year") +
    coord_cartesian(xlim = c(2018, 2026), ylim = c(0.004, 0.020)) +
    theme(legend.position = 'none')


  ggsave(sf_prevalence, filename = here("output/09-10-18/plots/sf_prevalence.png"), dpi = 300, width = unit(7, 'in'), height = unit(5, 'in'))
  ggsave(ba_prevalence, filename = here("output/09-10-18/plots/ba_prevalence.png"), dpi = 300, width = unit(7, 'in'), height = unit(5, 'in'))


  # number of tests
  tests_df <- select(df, site, interv, years_out, tests) %>%
    group_by(site, years_out, interv) %>%
    summarize(
      cummean = mean(tests, na.rm = T),
      ci_high = quantile(tests, .975, na.rm = T),
      ci_low = quantile(tests, .025, na.rm = T)
    )

  tests_df <- tests_df %>%
    group_by(site, interv) %>%
    mutate(mean = c(cummean[[1]], cummean[2:length(cummean)] - cummean[1:(length(cummean)-1)]))

  sf_tests <-
    ggplot(
      filter(tests_df, site == 'SF'),
      aes(x = years_out + 2016,
          y = mean,
          # ymax = ci_high,
          # ymin = ci_low,
          fill = interv,
          color = interv)) +
    geom_line() +
    geom_text_repel(
      data =
        filter(tests_df, site == 'SF', years_out == 5),
      mapping =
        aes(x = 2016 + years_out + .05,  label = interv),
      size = 3,
      nudge_x = .7,
      direction = "y",
      segment.size = .2,
      hjust = 0
    ) +
    ggtitle("Screening Tests in San Francisco", subtitle = "Simulating 5 years of hypothetical interventions") +
    ylab("Number of Tests per Year in the Intervention Time Period") +
    xlab("Year") +
    coord_cartesian(xlim = c(2017, 2025)) +
    theme(legend.position = 'none')

  ba_tests <-
    ggplot(
      filter(tests_df, site == 'BA'),
      aes(x = years_out + 2017,
          y = mean,
          # ymax = ci_high,
          # ymin = ci_low,
          fill = interv,
          color = interv)) +
    geom_line() +
    geom_text_repel(
      data =
        filter(tests_df, site == 'BA', years_out == 5),
      mapping =
        aes(x = 2017 + years_out + .05,  label = interv),
      size = 3,
      nudge_x = .7,
      direction = "y",
      segment.size = .2,
      hjust = 0
    ) +
    ggtitle("Screening Tests in Baltimore", subtitle = "Simulating 5 years of hypothetical interventions") +
    ylab("Number of Tests per Year in the Intervention Time Period") +
    xlab("Year") +
    coord_cartesian(xlim = c(2018, 2026)) +
    theme(legend.position = 'none')


  ggsave(sf_tests, filename = here("output/09-10-18/plots/sf_tests.png"), dpi = 300, width = unit(7, 'in'), height = unit(5, 'in'))
  ggsave(ba_tests, filename = here("output/09-10-18/plots/ba_tests.png"), dpi = 300, width = unit(7, 'in'), height = unit(5, 'in'))

  # Cases Averted

  cases_averted <- filter(df, interv != "base_case", years_out == 5) %>%
    select(sim, site, interv, years_out, cases_averted) %>%
    group_by(site, interv, years_out) %>%
    summarize(
      mean = mean(cases_averted, na.rm=T),
      ci_high = quantile(cases_averted, .975, na.rm=T),
      ci_low = quantile(cases_averted, .0025, na.rm = T),
      uqr = quantile(cases_averted, 0.75, na.rm = T),
      lqr = quantile(cases_averted, 0.25, na.rm = T)
    )


  sf_cases_averted <- ggplot(
    data = filter(cases_averted, site == 'SF'),
    mapping = aes(x = interv,
                  middle = mean,
                  ymax = ci_high,
                  ymin = ci_low,
                  upper = uqr,
                  lower = lqr,
                  fill = interv
    )
  ) +
    geom_boxplot(stat = 'identity') +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 90)) +
    ggtitle("Number of Incident Cases Averted in San Francisco", subtitle = "Simulating 5 years of hypothetical interventions") +
    xlab("Intervention") +
    ylab("Cases Averted Compared to Base Case")


  ba_cases_averted <- ggplot(
    data = filter(cases_averted, site == 'BA'),
    mapping = aes(x = interv,
                  middle = mean,
                  ymax = ci_high,
                  ymin = ci_low,
                  upper = uqr,
                  lower = lqr,
                  fill = interv
    )
  ) +
    geom_boxplot(stat = 'identity') +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 90)) +
    ggtitle("Number of Incident Cases Averted in Baltimore", subtitle = "Simulating 5 years of hypothetical interventions") +
    xlab("Intervention") +
    ylab("Cases Averted Compared to Base Case")


  ggsave(sf_cases_averted, filename = here("output/09-10-18/plots/sf_cases_averted.png"), dpi = 300, width = unit(10, 'in'), height = unit(7, 'in'))
  ggsave(ba_cases_averted, filename = here("output/09-10-18/plots/ba_cases_averted.png"), dpi = 300, width = unit(10, 'in'), height = unit(7, 'in'))


  # number needed to test
  nnt <- filter(df, interv != "base_case", ! grepl("MSM", interv), years_out == 5) %>%
    select(sim, site, interv, years_out, nnt) %>%
    group_by(site, interv, years_out) %>%
    summarize(
      mean = mean(nnt, na.rm=T),
      ci_high = quantile(nnt, .975, na.rm=T),
      ci_low = quantile(nnt, .0025, na.rm = T),
      uqr = quantile(nnt, 0.75, na.rm = T),
      lqr = quantile(nnt, 0.25, na.rm = T)
    )


  sf_nnt <- ggplot(
    data = filter(nnt, site == 'SF'),
    mapping = aes(x = interv,
                  y = mean
    )
  ) +
    # geom_boxplot(stat = 'identity') +
    geom_point() +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 90)) +
    ggtitle("Number of Tests Needed to Avert One Case in San Francisco", subtitle = "Simulating 5 years of hypothetical interventions") +
    xlab("Intervention") +
    ylab("Number of Tests Needed to Avert One Case Compared to Base Case")


  ba_cases_averted <- ggplot(
    data = filter(cases_averted, site == 'BA'),
    mapping = aes(x = interv,
                  middle = mean,
                  ymax = ci_high,
                  ymin = ci_low,
                  upper = uqr,
                  lower = lqr,
                  fill = interv
    )
  ) +
    geom_boxplot(stat = 'identity') +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 90)) +
    ggtitle("Number of Incident Cases Averted in Baltimore", subtitle = "Simulating 5 years of hypothetical interventions") +
    xlab("Intervention") +
    ylab("Cases Averted Compared to Base Case")




}


construct_interventions_df <- function() {
  require(here)
  results_df_rows <- list()
  for (site in c("SF", "BA")) {
    results <- switch(
      site,
      BA = readRDS(here("output/09-10-18/simulate_interventions/BA_complete.rds"))[["BA"]],
      SF = readRDS(here("output/09-10-18/simulate_interventions/SF_complete.rds"))[["SF"]]
    )

    intervention_names <- names(results[[1]])

    for (i in 1:length(results)) {
      for (interv in intervention_names) {
        row_entry <- data.frame(
          sim = i,
          site = site,
          interv = interv,
          years_out = 1:length(results[[i]][[interv]][['prevalence']]),
          prevalence = results[[i]][[interv]][['prevalence']],
          tests = results[[i]][[interv]][['tests']],
          cases_averted = ifelse(
            is.null(results[[i]][[interv]][['cases_averted']]),
            NA,
            results[[i]][[interv]][['cases_averted']]
          ),
          nnt = ifelse(
            is.null(results[[i]][[interv]][['nnt']]),
            NA,
            results[[i]][[interv]][['nnt']]
          )
        )
        results_df_rows[[length(results_df_rows)+1]] <- row_entry
      }
    }
  }

  results_df <- do.call(rbind.data.frame, results_df_rows)
  saveRDS(results, here("output/09-10-18/interventions_df.rds"))
}


