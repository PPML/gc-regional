#' Run Calibration, Render Prior/Posterior Plots, and Save Outputs
calibrate_gcRegional <-
  function(site,
           niterations,
           output_directory,
           calibration_filename,
           plot_filename,
           sd.factor = 100) {

    # set up variables
    load_start(site)
    sd.prior <- gc_env$sd.prior
    theta <- gc_env$theta


    # run calibration
    calibration <- my_mcmc(
      target = dLogPosterior,
      init.theta = theta,
      proposal.sd = sd.prior / sd.factor,
      n.iterations = niterations
    )

    # default calibration_filename
    if (missing(calibration_filename))
      calibration_filename <- paste0(Sys.Date(), " calibration_", site, ".rds")
    if (missing(plot_filename))
      plot_filename <- paste0(Sys.Date(), " plots_", site, ".pdf")

    # include gc_testing in the saved calibration
    if (exists('gc_testing', envir = .GlobalEnv)) {
      calibration$gc_testing <- gc_testing
    }

    # save calibration
    saveRDS(calibration, file.path(output_directory, calibration_filename))

    # plot outputs
    plot_outputs_from_calibration_file(
      calibration_rds_file = file.path(output_directory, calibration_filename)
      ,site = site
      ,output_directory = output_directory
      ,output_file = plot_filename
    )

    return(invisible(list(
      file.path(output_directory, calibration_filename),
      file.path(output_directory, plot_filename)
    )))
  }
