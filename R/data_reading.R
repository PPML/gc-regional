require(here)
require(readxl)
require(dplyr)
require(magrittr)

#' Get Excel Sheets as a named list

list_excel_sheets <- function(path) {
  sheet_names <- readxl::excel_sheets(path)
  l <- lapply(sheet_names, function(x) {
    readxl::read_excel(path, sheet = x)
  })
  names(l) <- sheet_names
  return(l)
}

#' Read in the SSuN 2016 Data Set from Baltimore
#'
#' This function is called immediately after its construction
#' so that the BA16_data is available for the gcRegional package.
load_BA16_data <- function() {

  BA16_filepath <- system.file("extdata/Baltimore_GC_data.xlsx", package = 'gcRegional')
  BA16_data <- list_excel_sheets(BA16_filepath)

  attach(BA16_data)
  `SSuN, Cases in MSM` %<>% .[1:14, ]
  `SSuN,Site of infection` %<>% .[,1:8]
  `SSuN, Time to medical care` %<>% .[1:6,]
  `Case notifcations` %<>% .[,1:6]
  `Case notifcations, race` %<>% .[,1:6]
  `Interview completion rates` %<>% .[,1:6]
  detach(BA16_data)

  BA16_data$`cal targets` <- NULL
  BA16_data$`Case notifcations (2)` <- NULL
  BA16_data$`Case notifcations, race (2)` <- NULL

  assign(x = 'BA16_data', value = BA16_data, envir = .GlobalEnv)
}


#' Read in the SSuN 2017 Data Set from Baltimore
#'
#' This function is called immediately after its construction
#' so that the BA16_data is available for the gcRegional package.

load_BA17_data <- function() {
  BA17_filepath <- system.file("extdata/regional_models_data_request_180509_baltimore.xlsx", package = 'gcRegional')
  BA17_data <- list_excel_sheets(BA17_filepath)

  attach(BA17_data)
  `Case notifcations` %<>% .[,1:7]
  `Case notifcations, race` %<>% .[,1:7]
  `Interview completion rates` %<>% .[1:32,]
  detach(BA17_data)

  assign(x = 'BA17_data', value = BA17_data, envir = .GlobalEnv)
  observed_proportion_msm_cases <- BA17_data$`SSuN, Cases in MSM`$`MSM GC cases` + BA17_data$`SSuN, Cases in MSM`$`Non-MSM GC cases`
  assign(x = 'observed_proportion_msm_cases', value = observed_proportion_msm_cases, envir = gc_env)
}


load_SF16_data <- function() {
  SF16_filepath <- system.file("extdata/SF_data.xlsx", package = 'gcRegional')
  SF16_data <- list_excel_sheets(SF16_filepath)
  assign(x = 'SF16_data', value = SF16_data, envir = .GlobalEnv)
    observed_cases_msm <- as.integer(unlist(SF16_data$`cases in MSM`[4:17, 3]))
    observed_cases_not_msm <- as.integer(unlist(SF16_data$`cases in MSM`[4:17, 3]))
    observed_proportion_msm_cases <- observed_cases_msm + observed_cases_not_msm
    assign(x = 'observed_proportion_msm_cases', observed_proportion_msm_cases, envir = gc_env)
}
