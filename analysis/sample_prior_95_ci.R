# Sample Prior 95% confidence Interval and Plot

site <- 'SF'

gc_testing <- c(
  'increase_risk_all_pops'
  # 'population_denominators',
  # 'logistic_screening',
  # 'only_increasing_screening',
  # 'national_posterior_medians'
  # 'population_risk_mixing'
)

load_start(site)

gc_env$scalars <- c(nhanes_prev = 1/100, case_rates = 1/10)

trace <- t(sapply(seq(0.025, 0.975, 0.01), function(x) qprior_transf(x)))
colnames(trace) <- names(gc_env$theta)

if (check_variation("national_posterior_medians"))
  trace <- t(apply(trace, 1, replace_national_parameter_medians))

trace.burn.thin <- trace.burn <- trace

post.sample <- model_fits(trace, use_trace_without_sampling = T)
pred <- as.data.frame(post.sample$outputs)

gc_assign(pred)

plot_posteriors(output_directory = '~/Desktop/',
                filename = paste0(site, "_sample",  ".pdf"))
rm(gc_env)
