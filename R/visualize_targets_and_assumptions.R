#' Visualize Calibration Target Data
visualize_targets_and_assumptions <- function(file = paste0(tempfile(), ".pdf")) {
  require(gcRegional)
  require(ggplot2)
  stopifnot('site' %in% gc_env)

  # NHANES ----
  nhanes_population_sizes <-
    ggplot(
      data = gc_env$nhanes.updated.dat,
      mapping = aes(x = race, y = prev_denom, fill = race)) +
    geom_bar(stat = "identity") +
    facet_wrap(~sex+age) +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle("NHANES Population Sizes for Prevalence Denominators 99-08")

  nhanes_prevalent_cases <-
    ggplot(
      data = gc_env$nhanes.updated.dat,
      mapping = aes(x = race, y = prev_numer, fill = race)) +
    geom_bar(stat = "identity") +
    facet_wrap(~sex+age) +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle("NHANES Prevalent Cases 99-08")

  pdf(file = file, paper="USr")


  return(file)
}
