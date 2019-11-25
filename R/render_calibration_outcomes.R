#' Load Calibration RDS files
#'
#' Given a path, this function reads in each calibration file
#' and returns a list containing the data from each RDS
#' file contained in that path.
#' @param path A directory containing RDS files which are the output of my_mcmc
load_calibration_rds_files <- function(path) {
  require(coda)
  stopifnot(dir.exists(path))
  files <- list.files(path, full.names = T)
  files_shortnames <- list.files(path, full.names = F)
  rds_files <- files[grepl(pattern = ".rds", x = files_shortnames, ignore.case = T)]
  stopifnot(length(rds_files) > 0)
  lapply(rds_files, function(x) readRDS(x))
}

#' Extract Traces from Calibrations
#' @param calibrations A list of calibration outputs from my_mcmc
extract_traces <- function(calibrations) {
  stopifnot(class(calibrations) == 'list')
  stopifnot(length(calibrations) > 0)
  stopifnot(all(sapply(calibrations, function(x) 'trace' %in% names(x))))
  mcmc.list(lapply(calibrations, function(x) as.mcmc(x$trace)))
}

#' Burn and Thin the Traces from a Directory Containing Calibration Data Files
burn_and_thin_traces_from_path <- function(path, burn, thin) {
  require(fitR)
  traces <- extract_traces(load_calibration_rds_files(path))
  trace <- burnAndThin(trace = traces, burn = burn, thin = thin)
  do.call(rbind, trace)
}

#' Render Predictions and Plots from Calibration Files
predictions_and_plots_from_files <- function(ba_path, sf_path, output_path, burn, thin, sample_n) {
  for (site in c("BA", "SF")) {
    load_start(site)
    trace <- burn_and_thin_traces_from_path(switch(site, BA = ba_path, SF = sf_path), burn = burn, thin = thin)
    trace <- trace[sample.int(n = nrow(trace), size = sample_n, replace = F), ]
    gc_env$trace.burn.thin <- gc_env$trace.burn <- gc_env$trace <- trace
    post.sample <- model_fits(trace, use_trace_without_sampling = T)
    gc_env$post.sample <- post.sample
    saveRDS(post.sample, file.path(output_path, paste0(site, "_post_sample.rds")))
    pred <- as.data.frame(post.sample$outputs)
    gc_assign(pred)
    theta.list <- as.data.frame(post.sample$theta.list)
    plot_posteriors(output_directory = output_path, filename = paste0(site, "_calibration_plots.pdf"))
  }
}
