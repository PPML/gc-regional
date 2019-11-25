# We're interested in why the MSM interventions have negative 
# numbers needed to treat for semi-annual but not annual or 
# quarter annual interventions. 

load_interventions_data <- function() {
	readRDS(system.file("intervention_analysis", "interventions2.rds", package='gcRegional'))
}

filter_msm_interventions_data <- function(df) {
	df %>% filter(
		interv %in% c('base_case', 'MSM_1.5', 'MSM_2.5', 'MSM_3', 'MSM_3.5', 'MSM_4.5', 'MSM_5', 'MSM_annual', 'MSM_semi_annual',
									'MSM_quarter_annual'),
	population == 'All',
	(site == 'BA' & year >= 2017) | 
	(site == 'SF' & year >= 2016)) %>% 
	select(-c(population, population_size))
}

# We'll plot the following time trends to try and get a better 
# picture of what's going on: 
# - prevalent infections 
# - diagnosed cases 
# - cumulative infections
# - cases averted
# - additional tests

reformat_msm_interventions_data <- function() {
	# if (!exists('interventions')) { 
		interventions <- load_interventions_data()
	# }
  interventions_msm <- filter_msm_interventions_data(interventions)
	interventions_msm %<>% 
		mutate(infections = symptomatic_cases + asymptomatic_cases) %>% 
		select(-c(symptomatic_cases, asymptomatic_cases)) 
	
	incident_cases_interv_dcast <- interventions_msm %>% 
		select(site, sim, interv, year, incident_cases) %>% 
		dcast(site + sim  + year ~ interv)

	cases_averted_df <- 
		mutate(incident_cases_interv_dcast, 
			`MSM_1.5` = base_case - `MSM_1.5`,
			`MSM_2.5` = base_case - `MSM_2.5`,
			`MSM_3` = base_case - `MSM_3`,
			`MSM_3.5` = base_case - `MSM_3.5`,
			`MSM_4.5` = base_case - `MSM_4.5`,
			`MSM_5` = base_case - `MSM_5`,
			MSM_annual = base_case - MSM_annual,
			MSM_quarter_annual = base_case - MSM_quarter_annual,
			MSM_semi_annual = base_case - MSM_semi_annual,
			base_case = NA)

	cases_averted_melt <- 
		melt(cases_averted_df, 
				 id.vars = c('site', 'sim', 'year'),
				 variable.name = 'interv',
				 value.name = 'cases_averted')

	interventions_msm %<>% merge(cases_averted_melt)

	additional_tests_interv_dcast <- interventions_msm %>%
		select(site, sim, interv, year, number_screened) %>% 
		dcast(site + sim  + year ~ interv)

	additional_tests_df <- 
		mutate(additional_tests_interv_dcast, 
			`MSM_1.5` = `MSM_1.5` - base_case,
			`MSM_2.5` = `MSM_2.5` - base_case,
			`MSM_3` = `MSM_3` - base_case,
			`MSM_3.5` = `MSM_3.5` - base_case,
			`MSM_4.5` = `MSM_4.5` - base_case,
			`MSM_5` = `MSM_5` - base_case,
			MSM_annual = MSM_annual - base_case,
			MSM_quarter_annual = MSM_quarter_annual - base_case,
			MSM_semi_annual = MSM_semi_annual - base_case,
		  base_case = NA)	

	additional_tests_melt <- 
		melt(additional_tests_df, 
				 id.vars = c('site', 'sim', 'year'),
				 variable.name = 'interv',
				 value.name = 'additional_tests')

	interventions_msm %<>% merge(additional_tests_melt)

  saveRDS(interventions_msm, file.path(
    system.file(
      'intervention_analysis/',
      package = 'gcRegional'), 'interventions_msm.rds'))

  interventions_msm_mean <- interventions_msm %>%
		group_by(site, interv, year) %>% 
		summarize_if(is.numeric, mean, na.rm=T)

	interventions_msm_mean_melt <- interventions_msm_mean %>% 
		melt(id.vars = c('site', 'sim', 'interv', 'year'))

	interventions_msm_mean_melt$interv <-
		as.factor(interventions_msm_mean_melt$interv)

	levels(interventions_msm_mean_melt$interv) <- 
	c('base case', '1.5', '2.5', '3', '3.5', '4.5', '5', '1', '4', '2')

  interventions_msm_mean_melt$variable <- as.factor(interventions_msm_mean_melt$variable)
	levels(interventions_msm_mean_melt$variable) <- 
		c('diagnosed cases', 'cumulative incident cases', 'number of screening tests done',
			'prevalent infections', 'incident cases averted', 'number of additional tests')

  saveRDS(interventions_msm_mean_melt, file.path(
    system.file(
      'intervention_analysis', package='gcRegional'
    ), 'interventions_mean_melt.rds'))

}

plot_msm_interventions <- function() { 

  interventions_msm_mean_melt <- 
    readRDS(file.path(
      system.file(
        'intervention_analysis', package='gcRegional'
      ), 'interventions_mean_melt.rds'))

	sf_plot <- ggplot(filter(interventions_msm_mean_melt, site=='SF'),
				 aes(x=year,
						 y=value,
						 group=interv,
						 color=interv)) + 
		geom_line() + 
		facet_wrap(~variable, scales='free_y') + 
		labs(color = 'screening rate per year in MSM')


	ba_plot <- ggplot(filter(interventions_msm_mean_melt, site=='BA'),
				 aes(x=year,
						 y=value,
						 group=interv,
						 color=interv)) + 
		geom_line() + 
		facet_wrap(~variable, scales='free_y') + 
		labs(color = 'screening rate per year in MSM')

  return(list(ba = ba_plot, sf = sf_plot))

	# interventions_msm_ci_high <- interventions_msm %>%
	# 	group_by(site, interv, year) %>% 
	# 	summarize_if(is.numeric, quantile, 0.975, na.rm=T)
			
	# interventions_msm_ci_low <- interventions_msm %>%
	# 	group_by(site, interv, year) %>% 
	# 	summarize_if(is.numeric, quantile, 0.025, na.rm=T)

	# interventions_msm_ci <- 
	# 	merge(
	# 		interventions_msm_ci_high,
	# 		interventions_msm_ci_low, 
	# 		by=c('site', 'interv', 'year'),
	# 		suffixes = c('.ci_high', '.ci_low'))

	# interventions_msm_ci_melt <- interventions_msm_ci %>% 
	# 	select(-c('sim.ci_high', 'sim.ci_low')) %>% 
	# 	melt(id.vars = c('site', 'interv', 'year')) %>% 
	# 	tidyr::separate(col='variable', into=c('variable', 'stat'), 
	# 									sep='\\.', fill='left') %>%
	# 	dcast(site + interv + year + variable ~ stat)

	# ggplot(filter(interventions_msm_ci_melt, variable %in% c('additional_tests', 'cases_averted'), site=='BA'),
	# 			 aes(x=year,
	# 					 ymax=ci_high,
	# 					 ymin=ci_low,
	# 					 group=interv,
	# 					 fill=interv)) + 
	# 	geom_ribbon(alpha=0.25) + 
    # # geom_line(aes(y=ci_low), linetype='dashed') + 
    # # geom_line(aes(y=ci_high), linetype='dashed') + 
	# 	geom_line(
	# 		data = filter(interventions_msm_mean_melt,
	# 									variable %in% c('additional_tests', 'cases_averted'),
	# 									site == 'BA'),
	# 		mapping = aes(x=year, y=value, group=interv, color=interv,
	# 									ymin=NULL,ymax=NULL)) + 
	# 	facet_wrap(~variable, scales='free_y')

	# interventions_msm_melt <- interventions_msm %>% 
	# 	melt(id.vars = c('site', 'sim', 'interv', 'year'))

		# facet_grid(interv~variable)
	
	# %>% 
		# mutate(incidents_1_year_ago = lag(incident_cases, 4),
		# 				 screens_1_year_ago = lag(diagnosed_cases, 4),
		# 				 nnt = (number_screened - screens_1_year_ago) / 
		# 					 (incident_cases - incidents_1_year_ago)) %>% 
		# select(-c(incidents_1_year_ago, screens_1_year_ago))

}
