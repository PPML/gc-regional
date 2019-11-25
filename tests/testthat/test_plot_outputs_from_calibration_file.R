test_that("plot_outputs_from_calibration_file produces a file given proper input.", {
  output_dir <- tempdir()
  output_file <- paste0(basename(tempfile()), '.pdf')
  filepath <- system.file(package = 'gcRegional', 'data/2018-11-16 BA-bezier-calibration.rds')
  if (! dir.exists(output_dir)) dir.create(output_dir)

  suppressWarnings(
  plot_outputs_from_calibration_file(
    calibration_rds_file = filepath
    ,site = "BA"
    ,burnin = 1000
    ,thin = 10
    ,sample_size = 2
    ,output_directory = output_dir
    ,output_file = output_file
  )
  )

  expect_true(file.exists(file.path(output_dir, output_file)))
})
