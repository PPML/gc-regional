# Let's try and optimize the likelihood within a 99% confidence intervals of the prior

library(magrittr)
library(gcRegional)

N_loops <- 20


site <-
  if (on_the_cluster()) {
    commandArgs(trailingOnly = T)[[1]]
  } else "BA" # default to one as a test
# and use commandLineArgs() if we're running
# on the cluster.

# minimum dLogPosterior acceptable for starting an optimization chain:
min_acceptable <- -1000
nth_sim <- if (on_the_cluster()) {
  nth_sim <-as.numeric(commandArgs(trailingOnly=T)[[2]])
} else 1
set.seed(nth_sim)

# output_directory <-
#   if(on_the_cluster()) "~/2018/October/27/optim-w-medians/" else "~/Documents/gcRegional/output/09-23-18/optim/"

gc_testing <- c(
  'population_denominators',
  # 'logistic_screening',
  # 'only_increasing_screening',
  'national_posterior_medians'
  # 'population_risk_mixing'
)


rm(gc_env) # remove any data lurking from a previous simulation

load_start(site)

gc_env$scalars <- c(nhanes_prev = 1/100, case_rates = 1/10, p.symp = 1/10)


theta_list <- list()
optim_list <- list()


# sample_until_acceptable <- function() {
#   theta <- rprior_transf()
#   accepted = F
#   while(! accepted) {
#     theta <- rprior_transf()
#     d <- dLogPosterior(theta)
#     accepted = is.finite(d) && d > min_acceptable
#   }
#   return(theta)
# }
# theta <- sample_until_acceptable()

trace <-
  switch(site, SF = readRDS("~/Documents/gcRegional/output/10-19-18/SF_optim_trace_top_15.rds"),
         BA = readRDS("~/Documents/gcRegional/output/10-19-18/BA_optim_trace_top_15.rds"))

theta <- unlist(trace[1,])
if (check_variation('population_risk_mixing') && ! 'logit.epsilon' %in% names(theta))
  theta['logit.epsilon'] <- logit(rbeta(1,1.1,1.1))

if (check_variation("national_posterior_medians"))
  theta <- remove_national_parameters(theta)

lik <- function(theta) {
  tryCatch({
  pred <- prediction_epi(theta)
  lik <- population_weighted_likelihood(pred)
  lik$p.symp <- lik$p.symp * 84 / 5
  lik$p.msm <- lik$p.msm * 84 / 16
  ll <- sum(sapply(lik, sum))
  if (! is.finite(ll)) return(-1e32)
  print(ll); return(ll)
  },
  error = function(x) -1e32
  )
}

out <- optim(
  par = theta,
  fn = lik,
  method = "Nelder-Mead",
  control = list(fnscale = -1)
)

theta <- out$par

trace <- do.call(rbind, lapply(1:10, function(x) theta))

if (check_variation('national_posterior_medians'))
  trace <- t(apply(trace, 1, replace_national_parameter_medians))

trace.burn.thin <- trace.burn <- trace

post.sample <- model_fits(trace, use_trace_without_sampling = T)
pred <- as.data.frame(post.sample$outputs)

gc_assign(pred)

saveRDS(object = theta, paste0("~/Desktop/", site, "_theta_lik_optim.rds"))

plot_posteriors(output_directory = "~/Desktop/",
                filename = paste0(site, "_likelihood_optim",  ".pdf"))
