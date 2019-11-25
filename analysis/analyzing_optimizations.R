# Analyzing Optimizations
library(magrittr)
# library(gcRegional)
devtools::load_all(".")


output_dir <- "~/Documents/gcRegional/output/01-07-19/"
setwd(output_dir)

get_len <- function(f) { length(readRDS(f)) }
# get_val <- function(f) { a <- readRDS(f); a[[length(a)]][['value']] }
get_val <- function(f) {
  a <- readRDS(f)
  vals <- sapply(1:length(a), function(i) a[[i]][['value']])
  return(max(vals))
}

sf_ls <- grep("SF_optim_list", list.files("optim/", full.names=T), value=T)
ba_ls <- grep("BA_optim_list", list.files("optim/", full.names=T), value=T)

sf_lengths <- sapply(sf_ls, get_len)
ba_lengths <- sapply(ba_ls, get_len)
sf_vals <- sapply(sf_ls, get_val)
ba_vals <- sapply(ba_ls, get_val)

get_theta <- function(f) {
  a <- readRDS(f)
  vals <- sapply(1:length(a), function(i) a[[i]][['value']])
  max_val <- which(vals == max(vals))
  return(a[[max_val[[1]]]][['par']])
}

ba_trace <- t(sapply(ba_ls, get_theta))
sf_trace <- t(sapply(sf_ls, get_theta))

sf_vals_sort <- sort(sf_vals)
top_sf_vals <- sf_vals_sort[(length(sf_vals_sort)-5):length(sf_vals_sort)]
top_sf_vals_idx <- which(sf_vals %in% top_sf_vals)
ba_vals_sort <- sort(ba_vals)
top_ba_vals <- ba_vals_sort[(length(ba_vals_sort)-5):length(ba_vals_sort)]
top_ba_vals_idx <- which(ba_vals %in% top_ba_vals)


ba_trace_top <- ba_trace[top_ba_vals_idx,]
sf_trace_top <- sf_trace[top_sf_vals_idx,]

saveRDS(ba_trace_top, paste0(output_dir, "BA_optim_trace_top.rds"))
saveRDS(sf_trace_top, paste0(output_dir, "SF_optim_trace_top.rds"))

gc_testing <-
  c(
    'population_denominators',
    'increase_risk_all_pops',
    'msm_incidence_likelihood',
    'separate_increased_risk',
    'only_use_nhanes_women',
    'only_increase_hetero_high_activity_transmission'
  )
gc_env$scalars <- c(nhanes_prev=1/100, case_rates = 1/10)


for (site in c("SF", "BA")) {
  load_start(site)
  trace <- switch(site, SF = sf_trace_top, BA = ba_trace_top)
  if (check_variation('national_posterior_medians')) {
    trace <- t(apply(trace, 1, replace_national_parameter_medians))
  }
  # colnames(trace) <- names(gc_env$theta)
  trace.burn.thin <- trace.burn <- trace

  post.sample <- model_fits(trace, use_trace_without_sampling = T)
  pred <- as.data.frame(post.sample$outputs)

  gc_assign(pred)

  plot_posteriors(output_directory = output_dir,
                  filename = paste0(site, "_optim_trace_top",  ".pdf"))
  rm(gc_env)

}
