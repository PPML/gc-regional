devtools::load_all()
library(here)

gc_testing <- c(
  'population_denominators',
  'msm_incidence_likelihood',
  'increase_risk_all_pops',
  'separate_increased_risk',
  'only_use_nhanes_women'
)

for (site in c('SF', 'BA')) {
	load_start(site)
	post.sample <- readRDS(here('inst/calibration_outcomes', paste0(site, '_posterior_sample.rds')))
	pred <- as.data.frame(post.sample$outputs)
	trace <- trace.burn <- trace.burn.thin <- as.data.frame(post.sample$theta.list)
	plots <- list_target_and_fit_plots()
	layout <- layout_target_and_fit_plots(plots)
	layout
	pdf(file=here('inst/calibration_outcomes/', paste0(site, '_calibration_plots.pdf')), paper='USr', width=11, height=8.5)
	print(layout)
	dev.off()
	plot_posteriors(filename=paste0(site, '_posterior_samples_2.pdf'), output_directory=here('inst/calibration_outcomes/'))
}
