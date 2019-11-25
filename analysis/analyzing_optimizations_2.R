# Get Optima from All Optimization Simulations as a Data Frame

library(gcRegional)


gc_testing <- c(
  'population_denominators'
  # 'only_increasing_screening',
  # 'logistic_screening'
)

output_directory <-
  "~/Documents/gcRegional/output/09-24-18/"


for (site in c("SF", "BA")) {

  rm(gc_env)
  load_start(site)
  gc_env$scalars <- c(nhanes_prev = 1/100, case_rates = 1/10)

  f_list <-
    grep(
      paste0(site, "_optim_list"),
      # list.files(paste0(output_directory, "optim-bezier/"), full.names = T),
      list.files(output_directory, full.names = T),
      value = T)

  trace <- do.call(rbind, lapply(f_list, function(f) {
    optim_out <- readRDS(f)
    t(sapply(1:length(optim_out), function(i) {
      return(optim_out[[i]][['par']])
      }))
  }))

  colnames(trace) <- names(gc_env$theta)

  trace_max <- apply(trace, 2, max)
  trace_min <- apply(trace, 2, min)

  dLogPosterior_vals <- unlist(sapply(f_list, function(f) {
    optim_out <- readRDS(f)
    sapply(1:length(optim_out), function(i) optim_out[[i]][['value']])
  }), use.names = F)

  convergence_vals <- unlist(sapply(f_list, function(f) {
    optim_out <- readRDS(f)
    sapply(1:length(optim_out), function(i) optim_out[[i]][['convergence']])
  }), use.names = F)

  dLogmax <- which(dLogPosterior_vals == max(dLogPosterior_vals[as.logical(convergence_vals)]))
  df <-
    data.frame(
      variable = names(gc_env$theta),
      min = apply(trace[as.logical(convergence_vals),], 2, min),
      max = apply(trace[as.logical(convergence_vals),], 2, max),
      y = trace[dLogmax, ],
      y_high = trace[dLogmax, ] + as.numeric(apply(trace[as.logical(convergence_vals), ], 2, sd)),
      y_low = trace[dLogmax, ] - as.numeric(apply(trace[as.logical(convergence_vals), ], 2, sd))
    )



  # keeper_threshold <- quantile(dLogPosterior_vals, .5)

  # keepers <- which(dLogPosterior_vals != -1e32)
  keepers <- which(convergence_vals == 0)

  trace <- trace[keepers, ]
  dLogPosterior_vals <- dLogPosterior_vals[keepers]

  keepers <- sort(dLogPosterior_vals)[(length(dLogPosterior_vals)-5):length(dLogPosterior_vals)]
  keepers <- which(dLogPosterior_vals %in% keepers)
  trace <- trace[keepers, ]
  dLogPosterior_vals <- dLogPosterior_vals[keepers]



  saveRDS(trace, paste0(output_directory, site, "_optim_trace.rds"))
  saveRDS(dLogPosterior_vals, paste0(output_directory, site, "_optim_vals.rds"))

  trace.burn.thin <- trace.burn <- trace

  post.sample <- model_fits(trace, use_trace_without_sampling = T)
  pred <- as.data.frame(post.sample$outputs)

  gc_assign(pred)

  plot_posteriors(output_directory = output_directory,
                  filename = paste0(site, "_optim_trace_", ".pdf"))




}
