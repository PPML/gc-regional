
#' Population Indicator Matrix
#'
#' This matrix's column names contain the different demographic categories,
#' i.e. race-ethnicity, gender, activity rate, and age category.
#' Its rows' indices refer to compartments in the model.
#' For example, one can query the matrix for what demographics
#' a compartment represents by looking at popmatrix[i, ], or
#' similarly one can look at which compartments represent a
#' particular demographic by looking at popmatrix[, i].
#'
#' Be careful not to use the Female MSM! They are in the model
#' for consistency, so that we can loop through populations
#' using a selection format like (i*j*k*l), but they are
#' always at 0 population size!
#' @name popmatrix
gc_env$popmatrix <-
structure(list(Male = c("1", "1", "1", "1", "1", "1", "1", "1",
"1", "1", "1", "1", "1", "1", "1", "1", "", "", "", "", "", "",
"", "", "", "", "", "", "", "", "", ""), Female = c("", "", "",
"", "", "", "", "", "", "", "", "", "", "", "", "", "1", "1",
"1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1",
"1"), Black = c("1", "1", "1", "1", "", "", "", "", "", "", "",
"", "", "", "", "", "1", "1", "1", "1", "", "", "", "", "", "",
"", "", "", "", "", ""), Other = c("", "", "", "", "1", "1",
"1", "1", "", "", "", "", "", "", "", "", "", "", "", "", "1",
"1", "1", "1", "", "", "", "", "", "", "", ""), Hispanic = c("",
"", "", "", "", "", "", "", "1", "1", "1", "1", "", "", "", "",
"", "", "", "", "", "", "", "", "1", "1", "1", "1", "", "", "",
""), MSM = c("", "", "", "", "", "", "", "", "", "", "", "",
"1", "1", "1", "1", "", "", "", "", "", "", "", "", "", "", "",
"", "1", "1", "1", "1"), Young = c("1", "1", "", "", "1", "1",
"", "", "1", "1", "", "", "1", "1", "", "", "1", "1", "", "",
"1", "1", "", "", "1", "1", "", "", "1", "1", "", ""), Old = c("",
"", "1", "1", "", "", "1", "1", "", "", "1", "1", "", "", "1",
"1", "", "", "1", "1", "", "", "1", "1", "", "", "1", "1", "",
"", "1", "1"), `High Activity` = c("", "1", "", "1", "", "1",
"", "1", "", "1", "", "1", "", "1", "", "1", "", "1", "", "1",
"", "1", "", "1", "", "1", "", "1", "", "1", "", "1"), `Low Activity` = c("1",
"", "1", "", "1", "", "1", "", "1", "", "1", "", "1", "", "1",
"", "1", "", "1", "", "1", "", "1", "", "1", "", "1", "", "1",
"", "1", "")), .Names = c("Male", "Female", "Black", "Other",
"Hispanic", "MSM", "Young", "Old", "High Activity", "Low Activity"
), row.names = c(NA, -32L), class = "data.frame")



# We need the index_to_pop_names and pop_names_to_index functions
# to convert back and forth from something like "Black MSM" to a
# list of columns in a given compartment in the gcSim output data.
index_to_pop_names <- function(index) {
colnames(gc_env$popmatrix)[which(gc_env$popmatrix[index, ] == 1)]
}

pop_names_to_index <- function(names) {
which(apply(gc_env$popmatrix, 1, function(x) all(x[names] == 1)))
}


gc_env$popnames <- list(
i = c("Black", "Hispanic", "Other", "MSM"),
j = c("Male", "Female"),
k = c("Low Activity", "High Activity"),
l = c("Young", "Old")
)

#' Select Population Names by their Indices
#'
#' This function converts a subset of c('i', 'j', 'k', 'l')
#' into a list of all possible intersections between each
#' of the designated populations.
#'
#' @param ijkl_subset Choose a subset of c('i', 'j', 'k', 'l') to pick subpopulations split by
#' the characteristics corresponding to each of these letters. For the correspondence between
#' these letters and the characteristics, see the gc_env$popnames object.
#'
#' @example select_populations(c('i', 'j'))
#' #>         i      j
#' #> 1    Black   Male
#' #> 2 Hispanic   Male
#' #> 3    White   Male
#' #> 4      MSM   Male
#' #> 5    Black Female
#' #> 6 Hispanic Female
#' #> 7    White Female
select_populations <- function(ijkl_subset = c('i', 'j', 'k', 'l')) {
selected_populations <- lapply(ijkl_subset, function(x) gc_env$popnames[[x]])
# selected_pop <- list()
# for (i in 1:nrow(selected_populations)) {
# selected_pop[[length(selected_pop)+1]] <- selected_populations[i, ]
# }
# selected_populations <- selected_pop; rm(selected_pop)
  selected_populations <- expand.grid(selected_populations)
  # selected_populations <- expand.grid(as.list(as.data.frame(selected_populations)))
  female_msm <-
    which(
      apply(selected_populations, 1, function(x)
        is.element("Female", x)) &
        apply(selected_populations, 1, function(x)
          is.element("MSM", x))
    )
  if (length(female_msm)!=0) {
    selected_populations <- selected_populations[-female_msm, ]
  }
  return(selected_populations)
}

#' Calculate the Percentage of Burden for Subpopulations
#'
#' We start with sol, the output of gcSim which is stored in gc_env by the
#' gcRegional::prediction_epi. sol is a dataframe which contains the
#' population sizes of each compartment at each timestep. For this function, we
#' are interested only in infections, and we subset sol for the
#' symptomatic/asymptomatic/all infections based on the which_infections
#' parameter. Using the ijkl_subset this function uses the
#' select_populations to group subpopulations and calculate the number
#' ofinfections experienced by a given subpopulation. Then the number of
#' infections of each subpopulation is divided by the total number of infections
#' to calculate the relative burden of each subpopulation.
#'
#' @param sol The output of the C++ gcSim function, stored by prediction_epi as gc_env$sol
#' @inheritParams select_populations
#' @param which_infections Specify "all" to tabulate both symptomatic and asymptomatic infections,
#' "y" for symptomatic, and "z" for symptomatic only.
#' @param relative If relative is TRUE, then each subpopulation is divided by the total number of
#' infections so that the output values represent each subpopulations percentage of the total
#' infections. If relative is FALSE, raw numbers of infections are output.
#' @param years Specify the vector of years to output prevalence data from.
#' The usual options for the years parameter are
#' 2010:2016  for the data period,
#' 2002:2016  for the calibration period, or
#' 1942:2016  for the entire simulation period.
tidy_prevalence_by_demographic <- function(
  sol,
  ijkl_subset = c("i", "j", "k", "l"),
  which_infections = "all",
  relative = TRUE,
  years = 2010:2016) {

  # 1942 == gc_env$start.year - gc_env$cal.start
  # -1 because of cpp's 0-indexing and R's 1-indexing ->
  # If we wanted to select 1942, the 0th year, we would
  # want to select the first row of sol.
  years <- years - 1941
  sol <- sol[years, ]

  # Collect the relevant infection data from the sol dataframe
  if (which_infections == "y") {
    y.index <- gc_env$y.index[1:28] # use 1:28 to exclude Female MSM
    infections <- sol[, 1+y.index]
  } else if (which_infections == "z") {
    z.index <- gc_env$z.index[1:28] # use 1:28 to exclude Female MSM
    infections <- sol[, 1+z.index]
  } else if (which_infections == "all") {
    y.index <- gc_env$y.index[1:28]
    z.index <- gc_env$z.index[1:28]
    infections <- sol[, 1+y.index]
    infections <- infections + sol[, 1+z.index]
  } else {
    stop("The which_infections parameter must be one of c('y', 'z', 'all)")
  }

  selected_populations <- select_populations(ijkl_subset)
  selected_indices <- list()
  for (i in 1:nrow(selected_populations)) {
    selected_indices[[i]] <-
      pop_names_to_index(as.character(unlist(selected_populations[i,])))
  }


  subpop_infections <-
    lapply(selected_indices, function(subpop_col_indices) {
      # Make sure we remove any consideration of "Female MSM"
      female_msm <- which(subpop_col_indices %in% 28:32)
      if (length(female_msm)!=0) {
        subpop_col_indices <- subpop_col_indices[-female_msm]
      }
      if (length(subpop_col_indices)>1) {
        return(rowSums(infections[, subpop_col_indices]))
      } else if (length(subpop_col_indices)==1) {
        return(infections[, subpop_col_indices])
      } else return(NULL)
    })

  null_entries_of_subpop_infections <-
    which(sapply(subpop_infections, is.null))
  if (length(null_entries_of_subpop_infections)!=0) {
    subpop_infections <- subpop_infections[-null_entries_of_subpop_infections]
    selected_populations <- selected_populations[-null_entries_of_subpop_infections, ]
  }
  subpop_infections <- as.data.frame(subpop_infections)


  # Rename the columns to represent the subpopulations
  colnames(subpop_infections) <-
    apply(selected_populations, 1, function(x)
      paste(x, collapse = " "))

  # Now we make the infection rates relative by dividing by
  # total infections at each time step
  total_infections <- rowSums(subpop_infections)
  # for (i in nrow(subpop_infections)) {
  #   subpop_infections[i, ] <- subpop_infections[i, ] / total_infections[[i]]
  # }
  if (relative) {
    subpop_infections <- subpop_infections / total_infections
  }
  # subpop_infections <- cbind(Year = 1:(nrow(subpop_infections)), subpop_infections)
  subpop_infections <- as.matrix(subpop_infections)

  # reshape2::melt the data for better formatting
  subpop_infections <- reshape2::melt(subpop_infections, varnames = (c("year", "demographic")))
  return(subpop_infections)
}


#' tidy_prevalence_dfs_for_many_simulations applies tidy_prevalence_by_demographic to many
#' simulations
#'
#'  @inheritParams tidy_prevalence_by_demographic
#'  @param sample_size The number of simulations to use for generating tidy prevalence data.
#'  interv An intervention plan that changes how screening is implemented after 2016.
#'  @param trace_sample A dataframe with rows of theta-values sampled from the trace of
#'  an MCMC calibration of the gcRegional model.
#'
#'  This function outputs a list of dataframes, each of which is generated individually
#'  by tidy_prevalence_by_demographic.
tidy_prevalence_dfs_for_many_simulations <-
  function(trace_sample,
           ijkl_subset = c('i', 'j', 'k', 'l'),
           relative = T,
           years = 2010:2016,
           which_infections = "all",
           interv = "base_case") {

  # We'll use each theta from trace_sample to run the gcSim deterministic compartmental model
  sol_list <- list()
  for (i in 1:nrow(trace_sample)) {
    theta <- setNames(as.numeric(trace_sample[i, ]), colnames(trace_sample))
    gcRegional:::prediction_epi(theta, interv=interv)
    sol_list[[length(sol_list)+1]] <- gc_env$sol
  }

  # Process each simulation into a dataframe with
  burden_df_list <-
    lapply(sol_list, function(x)
      tidy_prevalence_by_demographic(
        sol = x,
        ijkl_subset = ijkl_subset,
        relative = relative,
        years = years,
        which_infections = which_infections
      ))
  for (i in 1:length(burden_df_list)) {
    burden_df_list[[i]] <- cbind(sim = i, burden_df_list[[i]])
  }
  return(burden_df_list)
}

average_tidy_prevalence_dfs <- function(prevalence_dfs_list) {
  avg_burden <- prevalence_dfs_list[[1]]
  for (j in 1:nrow(prevalence_dfs_list[[1]])) {
    avg_burden[j, 4] <-
      mean(sapply(1:length(prevalence_dfs_list), function(i)
        prevalence_dfs_list[[i]][j, 4]), na.rm = T)
  }
  return(avg_burden)
}

#' This function converts y/z/all -> "Symptomatic ", "Asymptomatic ", ""
#' for the purpose of titling relative burden plots.
infection_str <- function(which_infections) {
  switch(which_infections,
         y = "Symptomatic ",
         z = "Asymptomatic ",
         all = "")
}


# Remember that this function takes the do.call(rbind, ...) version of the
# output of tidy_prevalence_dfs_for_many_simulations.
plot_faceted_prevalence <-
  function(prevalence_dfs,
           which_infections = "all",
           site = gc_env$site,
           years = 2010:2016,
           relative = T,
           scales = 'fixed') {

    prevalence_dfs <-
      dplyr::mutate(prevalence_dfs,
                    demographic = stringr::str_wrap(demographic, width = 15))
    prevalence_dfs$year <- years
    which_infections <- infection_str(which_infections)
    g <-
      ggplot(prevalence_dfs, aes(x = year, y = 100 * value, group = sim)) +
      geom_line(color = 'gray') + facet_wrap( ~ demographic, scales = scales) +
      stat_summary(
        aes(x = year, y = value * 100, group = 1),
        geom = "line",
        fun.y = "median",
        color = "black",
        size = 1
      ) +
      ggtitle(
        paste0(
          "Breakdown of Calibrated Gonorrhea ",
          which_infections,
          "Infections by Subpopulation in ",
          site
        )
      ) +
      xlab("Year") +
      ylab(
        ifelse(
          relative,
          "Percentage (%) of Infected Population",
          "Number (#) of Infected Individuals"
        )
      ) +
      labs(fill = 'Subpopulation') +
      theme(axis.text.x = element_text(angle = 90))
    return(g)
  }


# Remember that this only takes output from the average_tidy_prevalence_dfs
plot_stacked_area_prevalence <-
  function(avg_prevalence,
           which_infections = "all",
           site = gc_env$site,
           legend.position = "right",
           relative = T,
           years = 2010:2016) {

    which_infections <- infection_str(which_infections)
    avg_prevalence$year <- years
    if (relative) avg_prevalence$value <- avg_prevalence$value * 100
    g <-
      ggplot(avg_prevalence, aes(x = year, y = value, fill = demographic)) +
      geom_area() + geom_line(position = 'stack',
                              size = .2,
                              color = 'black') +
      ggtitle(paste0(
        "Breakdown of Calibrated Gonorrhea ",
        which_infections,
        "Infections by Subpopulation in ",
        site
      )) +
      xlab("Year") +
      ylab(
        ifelse(
          relative,
          "Percentage (%) of Total Infected Population",
          "Number (#) of Infected Individuals"
        )
      ) +
      labs(fill='Subpopulation') +
      # geom_text(data=dd, aes(x=year, y=cum*100, label=demographic), size=2) +
      theme(
        axis.text.x = element_text(angle = 90),
        legend.position = legend.position,
        legend.text = element_text(size = 10)
      )
    return(g)
  }


