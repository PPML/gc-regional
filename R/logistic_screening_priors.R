adjust_prior_for_logistic_screening <- function() {
  stopifnot(logistic_screening_variation())

  new_rows <- data.frame(
    parameter = c(paste0(unique(gc_env$screening_names), '.growth'), 'rep.symp.growth'),
    distribution = 'gamma',
    param1 = 2,
    param2 = 3,
    mean = 2 / 3,
    sd = sqrt(2 * 3 ^ 2),
    sd.prior = sqrt(2 * 3 ^ 2),
    transformation = 'log',
    mean.transf = log(2 / 3),
    ucl.transf = NA,
    sd.transf.1 = sd(log(rgamma(100000, 2, 3))),
    sd.transf = NA,
    sd1 = NA
  )
  gc_env$priors <- rbind(gc_env$priors, new_rows)
  gc_env$sd.prior <- gc_env$priors$sd.transf.1

  gc_env$priors <-
  gc_env$priors[
    ! grepl(paste0(paste0("*screen.", gc_env$screening_names, ".b"),collapse = "|"), gc_env$priors$parameter),
    ]

  gc_env$priors <-
  gc_env$priors[
    ! grepl(paste0(paste0("*screen.", gc_env$screening_names, ".c"),collapse = "|"), gc_env$priors$parameter),
    ]

  gc_env$theta <-
    gc_env$theta[
      ! grepl(paste0(paste0("*screen.", gc_env$screening_names, ".b"),collapse = "|"), names(gc_env$theta))
    ]
  gc_env$theta <-
    gc_env$theta[
      ! grepl(paste0(paste0("*screen.", gc_env$screening_names, ".c"),collapse = "|"), names(gc_env$theta))
    ]

  gc_env$priors[
    grepl(paste0(paste0("*screen.", gc_env$screening_names, ".d"), collapse = "|"), gc_env$priors$parameter), 'param1'] <- 1
  gc_env$priors[
    grepl(paste0(paste0("*screen.", gc_env$screening_names, ".d"), collapse = "|"), gc_env$priors$parameter) , 'param2'] <- 2


}
