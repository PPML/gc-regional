
save_national_natural_history_posterior <- function() {
  nat_posterior_sample <- readRDS("~/Documents/gc-regional-project/data/2018-06-25 national trace burn thin sample.rds")
  natural_history_names <-
    c(
      "logit.b.f",
      "logit.b.m",
      "logit.b.msm",
      "log.dur.inf.symp.m",
      "log.dur.inf.symp.msm",
      "log.dur.inf.symp.f",
      "log.dur.inf.asymp.m",
      "log.dur.inf.asymp.f",
      "log.dur.inf.asymp.msm",
      "logit.symp.m",
      "logit.symp.f",
      "logit.symp.msm"
    )
  natural_history <- nat_posterior_sample[, natural_history_names]
  saveRDS(natural_history, system.file("data/national_natural_history_parameters.rds", package = 'gcRegional'))
}


fit_national_natural_history_parameters <- function() {
  require(dplyr)
  nat_hist <- readRDS(system.file("data/national_natural_history_parameters.rds", package = 'gcRegional'))
  load_priors()
  param_names <- names(nat_hist)
  param_names <- strsplit(param_names, "\\.")
  transforms <- sapply(param_names, `[[`, 1)
  param_names <- # get the transformation for each column
    sapply(param_names, function(x) paste(x[2:length(x)], collapse = "."))
  colnames(nat_hist) <- param_names # rename columns without transformations
  logit_cols <- which(transforms == 'logit')
  log_cols <- setdiff(1:length(transforms), logit_cols)
  nat_hist <- mutate_at(nat_hist, logit_cols, ilogit) # undo logit transforms
  nat_hist <- mutate_at(nat_hist, log_cols, exp) # undo log transforms

  fits <- lapply(1:length(param_names), function(x) {
    param_row <- which(gc_env$priors$parameter == param_names[[x]])
    fun <- as.character(gc_env$priors[[param_row, 2]])
    startvals <- switch(
      fun,
      gamma = list(shape = gc_env$priors[param_row, 3], rate = gc_env$priors[param_row, 4]),
      beta = list(shape1 = gc_env$priors[param_row, 3], shape2 = gc_env$priors[param_row, 4]),
      normal = NULL
    )
    fit <- MASS::fitdistr(
      x = nat_hist[, x],
      densfun = fun,
      start = startvals
    )
    return(list(
      fit = fit,
      param = param_names[[x]],
      fun = fun))
  })

  saveRDS(fits, file.path(system.file("data", package = 'gcRegional'), "national_natural_history_fits.rds"))
}
