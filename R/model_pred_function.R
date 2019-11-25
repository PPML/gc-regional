###############################################################################################
### function to update estimated model parameters, run the model, and print desired outputs ###
###############################################################################################

aaply <- plyr::aaply

prediction_epi <- function(theta, interv = "base_case") {
  e <- create_gcSim_param_env(theta, interv = interv)
  gc_env$e <- e
  out.cpp <- run_gcSim_with_environment(e)
  output_statistics <- get_simulation_statistics(e)
  return(output_statistics)
}

#' Create Parameters for Calling the gcSim cppFunction
create_gcSim_param_env <- function(theta, interv = "base_case") {
  # We're going to make an environment which will contain all the
  # variables necessary to run the gcSim function.
  e <- new.env()
  with(e, {
    m1 <- gc_env$m1
    m2 <- gc_env$m2
    m3 <- gc_env$m3
    m4 <- gc_env$m4
    f1 <- gc_env$f1
    f2 <- gc_env$f2
    f3 <- gc_env$f3
    tstep <- gc_env$tstep
    model.end <- gc_env$model.end
    intervention.model.end <- gc_env$intervention.model.end
    yinit <- gc_env$yinit
    n.sa <- gc_env$n.sa
    births <- gc_env$births
    births.sa <- gc_env$births.sa
    births.nsa <- gc_env$births.nsa
    aging <- gc_env$aging
    aging.nsa<- gc_env$aging.nsa
    cal.start <- gc_env$cal.start
    cal.end <- gc_env$cal.end
    cal_and_intervention_rows <- gc_env$cal_and_intervention_rows

    ### return log/logit transformed parameters back to untransformed values and put in required form for use in transmission model ###
    epsilon <- ilogit(c(theta["logit.epsilon.1"], theta["logit.epsilon.2"], theta["logit.epsilon.3"], theta["logit.epsilon.4"])) #mixing between high and low activity groups (black, other, Hispanic, MSM)
    if (check_variation('population_risk_mixing')) {
      # epsilon_pop is the population average
      epsilon_pop <- ilogit(theta['logit.epsilon'])
      # epsilon is a vector of deviations from the population average
      # We scale this to -1 to 1
      # then we scale the population rate for each demographic
      # using epsilon as a percentage of the distance to 0 or 1 depending
      # on if epsilon is non-negative or negative.

      # So if epsilon.1 is .5, we assume the epsilon_pop for demographic 1.
      # If epsilon.1 is .75, we translate that to 50%, and say that the
      # risk mixing for demographic 1 is the population average risk mixing
      # plus 50% of the distance to 100% risk mixing. So if epsilon_pop
      # were .2, and epsilon.1 is .5, we then set the risk mixing for
      # demographic 1 to .6.
      epsilon <- 2*(epsilon-.5)
      epsilon <- sapply(epsilon, function(x) {
        if (x >= 0) {
          epsilon_pop + (1-epsilon_pop)*x
        } else {
          epsilon_pop + epsilon_pop*x
        }
      })
    }

    c.min.param <- exp(c(theta["log.c.min.m.1"], theta["log.c.min.m.2"], theta["log.c.min.f.1"], theta["log.c.min.f.2"], theta["log.c.min.msm.1"], theta["log.c.min.msm.2"]))
    c.min <- c_min_fun(c.min.param) #produces a vector of minimum partner change rates

    ctrl.pt<-update_ctrl(theta)  #get internal control points of bezier curves
    screen <- construct_screening_matrix(theta, ctrl.pt) #construct screening matrix

    # Generate screening rates for the pre-calibration period equal to the start of the varying screening rates
    # and post-calibration period screening rates equal to the end of the varying screening rates.
    n_intervention_years <- max(model.end - nrow(screen) - (cal.start), 0)
    screen <-
      rbind(screen[rep(1, times = (cal.start)), ],
            screen,
            screen[rep(nrow(screen), times = n_intervention_years), ])

    # Update screening rates based on selected intervention
    intervention_rows <- (nrow(screen) - n_intervention_years):nrow(screen)

    create_screening_intervention <-
      function(
        screen_mat = screen,
        screening_intervention_rows = intervention_rows,
        population_names,
        high_activity_targeting = FALSE,
        new_rate,
        targeting_effectivity = c(high = .5, low = 0.1),
        additive = F,
        only_increase = T
      ) {
        # The last_data_period_row is needed to get the screening rate before interventions start.
        last_data_period_row <- min(screening_intervention_rows)-1

        # If the intervention has only_increase = T, then it will only increase screening rates.
        # That is, if the original screening rates are higher than the rate
        # specified by the intervention, the rate will remain unchanged, and if the rate
        # specified by the intervention is higher than the original rate, the screening rate will
        # be updated.
        new_rate_f <- if (only_increase) function(x, orig) max(x, orig) else function(x, orig) x

        if (high_activity_targeting) {
          # Get the low_activity and high_activity column indices and screening rates.
          low_activity_cols <- pop_names_to_index(c(population_names, "Low Activity"))
          high_activity_cols <- pop_names_to_index(c(population_names, "High Activity"))
          low_activity_screen_rate <- screen_mat[last_data_period_row, low_activity_cols]
          high_activity_screen_rate <- screen_mat[last_data_period_row, high_activity_cols]

          # targeting_effectivity[['low']] says how many of the low activity individuals
          # are "accidentally" included in a targeted screening program, and the
          # targeting_effectivity[['high']] variable describes how many of the high
          # activity individuals are successfully targeted in their treatment.
          if (additive) {
            low_activity_screen_rate <-
              (1 - targeting_effectivity[['low']]) * low_activity_screen_rate +
              targeting_effectivity[['low']] * (low_activity_screen_rate + new_rate)
            high_activity_screen_rate <-
              (1 - targeting_effectivity[['high']]) * high_activity_screen_rate +
              targeting_effectivity[['high']] * (high_activity_screen_rate + new_rate)
          } else {
            low_activity_screen_rate <-
              (1 - targeting_effectivity[['low']]) * low_activity_screen_rate +
              targeting_effectivity[['low']] * new_rate
            high_activity_screen_rate <-
              (1 - targeting_effectivity[['high']]) * high_activity_screen_rate +
              targeting_effectivity[['high']] * new_rate
          }


          # Update screening
          for (i in screening_intervention_rows) {
            screen_mat[i, low_activity_cols] <-
              sapply(screen_mat[i, low_activity_cols], function(x) new_rate_f(x, low_activity_screen_rate))
            screen_mat[i, high_activity_cols] <-
              sapply(screen_mat[i, high_activity_cols], function(x) new_rate_f(x, high_activity_screen_rate))
          }
        } else {
          for (i in screening_intervention_rows) {
            screen_mat[i, pop_names_to_index(population_names)] <-
              sapply(screen_mat[i, pop_names_to_index(population_names)], function(x) new_rate_f(x, new_rate))
          }
        }

        # Return screening matrix
        return(screen_mat)
      }

    remove_ltfu_and_delay <-
      function() {
        switch(interv,
               'remove_ltfu_10pct' = {
                 screen[intervention_rows, ] <-
                   1/((1/screen[intervention_rows, ] )*(.9))
                 },
               'remove_ltfu_20pct' = {
                 screen[intervention_rows, ] <-
                   1/((1/screen[intervention_rows, ] )*(.8))
                 },
               'remove_f_msm_10pct_and_msw_20pct_ltfu' = {
                   F_and_MSM <- c(pop_names_to_index("Female"), pop_names_to_index("MSM"))
                   MSW <- setdiff(pop_names_to_index("Male"), pop_names_to_index("MSM"))
                   screen[intervention_rows, F_and_MSM] <-
                     1/((1/screen[intervention_rows, F_and_MSM] )*(.9))
                   screen[intervention_rows, MSW] <-
                     1/((1/screen[intervention_rows, MSW] )*(.8))
                 })
        return(screen)
      }


    switch(interv,
           "base_case" = {},
           "no_screening" = { screen[(2010-1942+1):nrow(screen),] <- 0 },
           "universal_annual" = { screen <- create_screening_intervention(population_names = c(), new_rate = 1)},
           "universal_semi_annual" = { screen <- create_screening_intervention(population_names = c(), new_rate = 2)},
           "high_activity_young_annual" = { screen <- create_screening_intervention(population_names = c("Young"), high_activity_targeting = TRUE, new_rate = 1) },
           "high_activity_young_semi_annual" = { screen <- create_screening_intervention(population_names = c("Young"), high_activity_targeting = TRUE, new_rate = 2) },
           "high_activity_young_quarter_annual" = { screen <- create_screening_intervention(population_names = c("Young"), high_activity_targeting = TRUE, new_rate = 4) },
           "young_annual" = { screen <- create_screening_intervention(population_names = c("Young"), new_rate = 1) },
           "young_semi_annual" = { screen <- create_screening_intervention(population_names = c("Young"), new_rate = 2) },
           "young_quarter_annual" = { screen <- create_screening_intervention(population_names = c("Young"), new_rate = 4) },
           "high_activity_MSM_annual" = { screen <- create_screening_intervention(population_names = c("MSM"), high_activity_targeting = TRUE, new_rate = 1) },
           "high_activity_MSM_semi_annual" = { screen <- create_screening_intervention(population_names = c("MSM"), high_activity_targeting = TRUE, new_rate = 2) },
           "high_activity_MSM_quarter_annual" = { screen <- create_screening_intervention(population_names = c("MSM"), high_activity_targeting = TRUE, new_rate = 4) },
           "MSM_annual" = { screen <- create_screening_intervention(population_names = c("MSM"), new_rate = 1) },
					 "MSM_1.5" = { screen <- create_screening_intervention(population_names = c("MSM"), new_rate = 1.5) },
           "MSM_semi_annual" = { screen <- create_screening_intervention(population_names = c("MSM"), new_rate = 2) },
					 "MSM_2.5" = { screen <- create_screening_intervention(population_names = c("MSM"), new_rate = 2.5) },
					 "MSM_3" = { screen <- create_screening_intervention(population_names = c("MSM"), new_rate = 3) },
					 "MSM_3.5" = { screen <- create_screening_intervention(population_names = c("MSM"), new_rate = 3.5) },
           "MSM_quarter_annual" = { screen <- create_screening_intervention(population_names = c("MSM"), new_rate = 4) },
					 "MSM_4.5" = { screen <- create_screening_intervention(population_names = c("MSM"), new_rate = 4.5) },
					 "MSM_5" = { screen <- create_screening_intervention(population_names = c("MSM"), new_rate = 5) },
           "high_activity_annual" = {  screen <- create_screening_intervention(population_names = c(), high_activity_targeting = TRUE, new_rate = 1) },
           "high_activity_semi_annual" = { screen <- create_screening_intervention(population_names = c(), high_activity_targeting = TRUE, new_rate = 2) },
           "high_activity_quarter_annual" = { screen <- create_screening_intervention(population_names = c(), high_activity_targeting = TRUE, new_rate = 3) },
           "female_young_annual" = { screen <- create_screening_intervention(population_names = c("Female", "Young"), new_rate = 1) },
           "female_young_semi_annual" = { screen <- create_screening_intervention(population_names = c("Female", "Young"), new_rate = 2) },
           "female_young_quarter_annual" = { screen <- create_screening_intervention(population_names = c("Female", "Young"), new_rate = 4) },
           "high_activity_female_young_annual" = { screen <- create_screening_intervention(population_names = c("Female", "Young"), high_activity_targeting = TRUE, new_rate = 1) },
           "high_activity_female_young_semi_annual" = { screen <- create_screening_intervention(population_names = c("Female", "Young"), high_activity_targeting = TRUE, new_rate = 2) },
           "high_activity_female_young_quarter_annual" = { screen <- create_screening_intervention(population_names = c("Female", "Young"), high_activity_targeting = TRUE, new_rate = 4) },
           'remove_ltfu_10pct' = { screen <- remove_ltfu_and_delay() },
           'remove_ltfu_20pct' = { screen <- remove_ltfu_and_delay() },
           'remove_f_msm_10pct_and_msw_20pct_ltfu' = { screen <- remove_ltfu_and_delay() },
           '1year_high_activity_semi_annual' = {
             screen <-
               create_screening_intervention(
                 population_names = c(),
                 screening_intervention_rows = intervention_rows[[1]],
                 high_activity_targeting = TRUE,
                 new_rate = 2,
                 additive = T
               )
           },
           '1year_high_activity_semi_annual_lower_capture_rate' = {
             screen <-
               create_screening_intervention(
                 population_names = c(),
                 screening_intervention_rows = intervention_rows[[1]],
                 high_activity_targeting = TRUE,
                 new_rate = 2,
                 targeting_effectivity = c(high = .25, low = .12778),
                 additive = T
               )
           },
           '1year_high_activity_semi_annual_complete_capture_rate' = {
             screen <-
               create_screening_intervention(
                 population_names = c(),
                 screening_intervention_rows = intervention_rows[[1]],
                 high_activity_targeting = TRUE,
                 new_rate = 2,
                 targeting_effectivity = c(high = 1, low = .04444),
                 additive = T
               )
           },
           '1year_high_activity_semi_annual_0_10th_hr_pt1555_lr' = {
             screen <-
               create_screening_intervention(
                 population_names = c(),
                 screening_intervention_rows = intervention_rows[[1]],
                 high_activity_targeting = TRUE,
                 new_rate = 2,
                 targeting_effectivity = c(high = 0, low = .1555556),
                 additive = T
               )
           },
           '1year_high_activity_semi_annual_1_10th_hr_pt144_lr' = {
             screen <-
               create_screening_intervention(
                 population_names = c(),
                 screening_intervention_rows = intervention_rows[[1]],
                 high_activity_targeting = TRUE,
                 new_rate = 2,
                 targeting_effectivity = c(high = .1, low = .14444),
                 additive = T
               )
           },
           '1year_high_activity_semi_annual_2_10th_hr_pt133_lr' = {
             screen <-
               create_screening_intervention(
                 population_names = c(),
                 screening_intervention_rows = intervention_rows[[1]],
                 high_activity_targeting = TRUE,
                 new_rate = 2,
                 targeting_effectivity = c(high = .2, low = .13333),
                 additive = T
               )
           },
           '1year_high_activity_semi_annual_3_10th_hr_pt122_lr' = {
             screen <-
               create_screening_intervention(
                 population_names = c(),
                 screening_intervention_rows = intervention_rows[[1]],
                 high_activity_targeting = TRUE,
                 new_rate = 2,
                 targeting_effectivity = c(high = .3, low = .12222),
                 additive = T
               )
           },
           '1year_high_activity_semi_annual_4_10th_hr_pt111_lr' = {
             screen <-
               create_screening_intervention(
                 population_names = c(),
                 screening_intervention_rows = intervention_rows[[1]],
                 high_activity_targeting = TRUE,
                 new_rate = 2,
                 targeting_effectivity = c(high = .4, low = .11111),
                 additive = T
               )
           },
           '1year_high_activity_semi_annual_5_10th_hr_pt1_lr' = {
             screen <-
               create_screening_intervention(
                 population_names = c(),
                 screening_intervention_rows = intervention_rows[[1]],
                 high_activity_targeting = TRUE,
                 new_rate = 2,
                 targeting_effectivity = c(high = .5, low = .1),
                 additive = T
               )
           },
           '1year_high_activity_semi_annual_6_10th_hr_pt088_lr' = {
             screen <-
               create_screening_intervention(
                 population_names = c(),
                 screening_intervention_rows = intervention_rows[[1]],
                 high_activity_targeting = TRUE,
                 new_rate = 2,
                 targeting_effectivity = c(high = .6, low = .08889),
                 additive = T
               )
           },
           '1year_high_activity_semi_annual_8_10th_hr_pt066_lr' = {
             screen <-
               create_screening_intervention(
                 population_names = c(),
                 screening_intervention_rows = intervention_rows[[1]],
                 high_activity_targeting = TRUE,
                 new_rate = 2,
                 targeting_effectivity = c(high = .8, low = .066667),
                 additive = T
               )
           }
    )

    b <- ilogit(c(theta["logit.b.m"], theta["logit.b.f"], theta["logit.b.msm"])) #transmission rates
    gamma <- c(365/exp(theta["log.dur.inf.symp.m"]), 365/exp(theta["log.dur.inf.symp.msm"]), 365/exp(theta["log.dur.inf.symp.f"])) #duration of infectiousness if symptomatic (d)
    delta <- c(365/exp(theta["log.dur.inf.asymp.m"]), 365/exp(theta["log.dur.inf.asymp.msm"]), 365/exp(theta["log.dur.inf.asymp.f"])) #duration of infectiousness if asymptomatic (d)
    symp <- ilogit(c(theta["logit.symp.m"],theta["logit.symp.msm"],theta["logit.symp.f"])) #probability infection is symptomatic
    theta.param <- ilogit(c(theta["logit.theta.1"], theta["logit.theta.2"], theta["logit.theta.3"], theta["logit.theta.4"], theta["logit.theta.5"], theta["logit.theta.6"], theta["logit.theta.7"],logit(1))) #subpopulation assortativity
    #theta.param <- ilogit(c(theta["logit.theta.1"], theta["logit.theta.2"], theta["logit.theta.3"], logit(0.99), theta["logit.theta.5"], theta["logit.theta.6"], theta["logit.theta.7"],logit(1)))
    theta.param <- matrix(theta.param,ncol=2)  # col 1= males, col 2= females
    pi.all<- ilogit(c(rep(c(theta["logit.pi.m"],theta["logit.pi.f"]),3),theta["logit.pi.msm"])) #age assortativity
    rr.rep.symp.m <- ilogit(theta["logit.risk.rep.symp.m"]) #relative risk case is reported if sympomatic and male (other or Hispanic)
    rr.rep.symp.f <-ilogit(theta["logit.risk.rep.symp.f"]) #relative risk case is reported if sympomatic and female (other or Hispanic)
    rr.rep.symp.m1 <- ilogit(theta["logit.risk.rep.symp.m1"]) #relative risk case is reported if sympomatic and male (black)
    rr.rep.symp.f1 <-ilogit(theta["logit.risk.rep.symp.f1"]) #relative risk case is reported if sympomatic and female (black)

    rep.trend <- make_reporting_trend(theta, ctrl.pt)
    rep <-rbind(rep.trend[rep(1, times=(cal.start-1)),], rep.trend, rep.trend[rep((nrow(rep.trend)), times=max(model.end - nrow(rep.trend) - (cal.start - 1), 0)),]) #expand reporting prob to cover entire model run time, assume reporting pre-calibration period = reporting at start of cal
    rep.symp<-matrix(c(rep[,m1]*rr.rep.symp.m1,rep[,c(m2,m3,m4)]*rr.rep.symp.m,rep[,f1]*rr.rep.symp.f1, rep[,c(f2,f3,f3)]*rr.rep.symp.f), ncol=32) #reporting rate if symptomatic male or female
    rep.asymp<-rep[,1] #reporting prob is assumed to be invariant across age/sex/subpop when case is asymptomatic, so don't need separate estimates
    behav.trend<-behav_fun(ilogit(theta["logit.behav.lin"])) # behav time trend representing changing condom use/behaviour in MSM
    behav <-c(rep(0, times=(cal.start-1)), behav.trend, behav.trend[rep(length(behav.trend), times=max(model.end - length(behav.trend) - (cal.start - 1), 0))]) #expand to cover entire model run time, assume behav pre-calibration period = value at start of calibration

    if (exists('gc_testing', envir=.GlobalEnv) && 'DelayedMSMTransmissionIncrease' %in% gc_testing) {
      behav <- c(rep(1, times=(cal.start-1+8)), behav.trend)
    }
    if (exists('gc_testing', envir=.GlobalEnv) && 'zero_msm_transmission_increase' %in% gc_testing) {
      behav <- rep(1, times=(cal.start+cal.end))
    }

    if (check_variation("separate_increased_risk")) {
      behav.hetero <- behav_fun(ilogit(theta['logit.behav.hetero']))
      behav.hetero <-c(rep(0, times=(cal.start-1)), behav.hetero, behav.hetero[rep(length(behav.hetero), times=max(model.end - length(behav.hetero) - (cal.start - 1), 0))]) #expand to cover entire model run time, assume behav pre-calibration period = value at start of calibration
    } else {
      behav.hetero <- behav
    }


    rp.all <- exp(c(theta["log.rp.1.1.1.1"],theta["log.rp.1.1.2.1"],theta["log.rp.1.2.1.1"],theta["log.rp.1.2.2.1"],  #relative rates of partner change in different population groups
                    theta["log.rp.1.1.1.2"],theta["log.rp.1.1.2.2"],theta["log.rp.1.2.1.2"], theta["log.rp.1.2.2.2"],
                    theta["log.rp.2.1.1.1"],theta["log.rp.2.1.2.1"],theta["log.rp.2.2.1.1"],theta["log.rp.2.2.2.1"],
                    theta["log.rp.2.1.1.2"],theta["log.rp.2.1.2.2"],theta["log.rp.2.2.1.2"], theta["log.rp.2.2.2.2"],
                    log(1),theta["log.rp.3.1.2.1"],log(1),theta["log.rp.3.2.2.1"],
                    log(1),theta["log.rp.3.1.2.2"],log(1), theta["log.rp.3.2.2.2"],
                    log(1),theta["log.rp.4.1.2.1"],log(1),log(1),
                    log(1),theta["log.rp.4.1.2.2"],log(1), log(1)))
    pred <- mixing(epsilon, pi.all, theta.param, c.min, rp.all)  #calculate mixing matrix based on current parameter estimates
    part.all <- gc_env$part.all
    part.all.m <- gc_env$part.all.m
    part.all.f <- gc_env$part.all.f
    cm.list <- gc_env$cm.list
    age.dist.all<-prop.table(apply(part.all, c(3:5), sum), c(1,3)) #get proportion of partnerships with same age group, by sex
    age.dist.all<-c(diag(age.dist.all[,,1]),diag(age.dist.all[,,2]))
    age.dist.sub<-prop.table(apply(part.all, c(3:6), sum),c(1,3,4) )  #get proportion of partnerships with same age group, by subpopulation and sex
    d.m=c()
    d.f=c()
    for (i in 1:4){
      d.m[(2*i-1):(2*i)]<-diag(age.dist.sub[,,1,i])
      d.f[(2*i-1):(2*i)]<-diag(age.dist.sub[,,2,i])
    }
    age.dist.all <- c(age.dist.all, d.m[1:(3*2)], d.f[1:(3*2)],d.m[7:8])  #add in estimate for MSM
    s.dist.m <- part.all.m[1:3,1:3]/apply(part.all.m[1:3,1:3],1,sum) #proportion of partners of same or other subpopulation, M (excluding MSM)
    s.dist.f <- part.all.f[1:3,1:3]/apply(part.all.f[1:3,1:3],1,sum) #proportion of parters of same or other subpopulation, F
    pred.s.dist <- c(diag(s.dist.m),diag(s.dist.f)) #parterns of same subpopulation, M and F
    params <-list(b=b,symp=symp,gamma=gamma,delta=delta,rep.symp=rep.symp, rep.asymp=rep.asymp, screen=screen, behav=behav, behav.hetero=behav.hetero) #parameters used by transmission model
    quarterly <- F
    increase_risk_all_pops <- check_variation('increase_risk_all_pops')
  })

  return(e)
}


#' Run gcSim with Environment Variables
run_gcSim_with_environment <- function(e) {
  if (check_variation("only_increase_hetero_high_activity_transmission")) {
    with(e, {
    out.cpp <- gcSim2(
      params,
      tstep,
      cm.list,
      model.end,
      yinit,
      n.sa,
      births,
      births.sa,
      births.nsa,
      aging,
      aging.nsa,
      debug = FALSE,
      quarterly = quarterly,
      increase_risk_all_pops = increase_risk_all_pops
    )
    return(e)
    })
  } else {
    with(e, {
    out.cpp <- gcSim(
      params,
      tstep,
      cm.list,
      model.end,
      yinit,
      n.sa,
      births,
      births.sa,
      births.nsa,
      aging,
      aging.nsa,
      debug = FALSE,
      quarterly = quarterly,
      increase_risk_all_pops = increase_risk_all_pops
    )
    return(e)
    })
  }
}


#' Render Statistics from a Simulation
get_simulation_statistics <- function(e) {
  with(e, {
    sol <- out.cpp$out #save transmission model output
    gc_assign(sol)
    prev<-model_prev(sol)  #calculate desired outputs for calibration, function stored in distribution.param.est.R
    output <- list(pred.s.dist=pred.s.dist, age.dist.all=age.dist.all, prev=prev)
    return(output)
  })
}


construct_screening_matrix <- function(theta, ctrl.pt) {
  if (logistic_screening_variation()) {
    screening <-
      matrix(
        sapply(gc_env$screening_names, function(x) {
          c(ilogit(theta[paste0('logit.screen.', x, '.a')]),
            ilogit(theta[paste0('logit.screen.', x, '.d')]),
            exp(theta[paste0('log.', x, '.growth')])
        )
      }), ncol=3, byrow=T)
    t_vec <-seq(0,1, length.out=(gc_env$cal.period+1))
    screen_mat <- apply(screening, 1, function(x) logistic_growth(x[[1]], x[[2]], x[[3]], t_vec))
  } else if (exists('gc_testing', envir = .GlobalEnv) &&
      'flat_screening' %in% gc_testing) {
    screen.bezier <-
      matrix(c(
        rep(ilogit(theta["logit.screen.m1.1.a"]), 4), # male black young
        rep(ilogit(theta["logit.screen.m2.1.a"]), 4), # male black old
        rep(ilogit(theta["logit.screen.m1.a"]), 4), # male other young
        rep(ilogit(theta["logit.screen.m2.a"]), 4), # male other old
        rep(ilogit(theta["logit.screen.m1.a"]), 4), # male hispanic young
        rep(ilogit(theta["logit.screen.m2.a"]), 4), # male hispanic old
        rep(ilogit(theta['logit.screen.msm1.a']), 4), # male msm young
        rep(ilogit(theta['logit.screen.msm2.a']), 4), # male msm old
        rep(ilogit(theta["logit.screen.f1.1.a"]), 4), # female black young
        rep(ilogit(theta["logit.screen.f2.1.a"]), 4), # female black old
        rep(ilogit(theta["logit.screen.f1.a"]), 4), # female other young
        rep(ilogit(theta["logit.screen.f2.a"]), 4), # female other old
        rep(ilogit(theta["logit.screen.f1.a"]), 4), # female hispanic young
        rep(ilogit(theta["logit.screen.f2.a"]), 4) # female hispanic old
      ), ncol=4, byrow=TRUE)
    screen_mat <- apply(screen.bezier,1, bezier_fun) #get base annual screening rates over the calibration period
  } else if (exists('gc_testing', envir = .GlobalEnv) &&
             'linear_screening' %in% gc_testing) {

    t_vec <-seq(0,1, length.out=(gc_env$cal.period+1))

    screening <- matrix(c(
    ilogit(theta["logit.screen.m1.1.a"]),
    ilogit(theta["logit.screen.m1.1.d"]),
    ilogit(theta["logit.screen.m2.1.a"]),
    ilogit(theta["logit.screen.m2.1.d"]),
    ilogit(theta["logit.screen.m1.a"]),
    ilogit(theta["logit.screen.m1.d"]),
    ilogit(theta["logit.screen.m2.a"]),
    ilogit(theta["logit.screen.m2.d"]),
    ilogit(theta["logit.screen.m1.a"]),
    ilogit(theta["logit.screen.m1.d"]),
    ilogit(theta["logit.screen.m2.a"]),
    ilogit(theta["logit.screen.m2.d"]),
    ilogit(theta['logit.screen.msm1.a']),
    ilogit(theta['logit.screen.msm1.d']),
    ilogit(theta['logit.screen.msm2.a']),
    ilogit(theta['logit.screen.msm2.d']),
    ilogit(theta["logit.screen.f1.1.a"]),
    ilogit(theta["logit.screen.f1.1.d"]),
    ilogit(theta["logit.screen.f2.1.a"]),
    ilogit(theta["logit.screen.f2.1.d"]),
    ilogit(theta["logit.screen.f1.a"]),
    ilogit(theta["logit.screen.f1.d"]),
    ilogit(theta["logit.screen.f2.a"]),
    ilogit(theta["logit.screen.f2.d"]),
    ilogit(theta["logit.screen.f1.a"]),
    ilogit(theta["logit.screen.f1.d"]),
    ilogit(theta["logit.screen.f2.a"]),
    ilogit(theta["logit.screen.f2.d"])
    ), ncol=2, byrow=T)

    screen_mat <- apply(screening, 1, function(x) x[[1]]*(1-t_vec)+x[[2]]*t_vec)

  } else {
    bez.bc <-ctrl.pt[[1]] #internal control points for screening

    screen.bezier <-
      matrix(
        c(c(
          ilogit(c(theta["logit.screen.m1.1.a"],theta["logit.screen.m1.1.d"])),
          bez.bc["m1.1","b"],
          bez.bc["m1.1","c"],
          ilogit(c(theta["logit.screen.m2.1.a"],theta["logit.screen.m2.1.d"])),
          bez.bc["m2.1","b"],
          bez.bc["m2.1","c"]),
          rep(c(ilogit(c(theta["logit.screen.m1.a"],theta["logit.screen.m1.d"])),
                bez.bc["m1","b"],
                bez.bc["m1","c"],
                ilogit(c(theta["logit.screen.m2.a"],theta["logit.screen.m2.d"])),
                bez.bc["m2","b"],
                bez.bc["m2","c"]), 2),
          c(
            if (exists('gc_testing', envir = .GlobalEnv) &&
                'flat_msm_screening' %in% gc_testing) {
              # if 'flat_msm_screening' is turned on, just use the first parameter as all
              # bezier control points
              rep(ilogit(theta['logit.screen.msm1.a']), 4)
            } else {
              c(ilogit(c(theta["logit.screen.msm1.a"], theta["logit.screen.msm1.d"])),
                bez.bc["msm1", "b"],
                bez.bc["msm1", "c"])
            },
            if (exists('gc_testing', envir = .GlobalEnv) &&
                'flat_msm_screening' %in% gc_testing) {
              # if 'flat_msm_screening' is turned on, just use the first parameter as all
              # bezier control points
              rep(ilogit(theta['logit.screen.msm2.a']), 4)
            } else {
              c(ilogit(c(theta["logit.screen.msm2.a"], theta["logit.screen.msm2.d"])),
                bez.bc["msm2", "b"],
                bez.bc["msm2", "c"])
            }
          ),
          c(ilogit(c(theta["logit.screen.f1.1.a"],theta["logit.screen.f1.1.d"])),
            bez.bc["f1.1","b"],
            bez.bc["f1.1","c"],
            ilogit(c(theta["logit.screen.f2.1.a"],theta["logit.screen.f2.1.d"])),
            bez.bc["f2.1","b"],
            bez.bc["f2.1","c"]),
          rep(c(ilogit(c(theta["logit.screen.f1.a"],theta["logit.screen.f1.d"])),
                bez.bc["f1","b"],
                bez.bc["f1","c"],
                ilogit(c(theta["logit.screen.f2.a"],theta["logit.screen.f2.d"])),
                bez.bc["f2","b"],
                bez.bc["f2","c"]),2)),
        ncol=4, byrow=TRUE)
    screen_mat <- apply(screen.bezier,1, bezier_fun) #get base annual screening rates over the calibration period

  }

  rr.screen.s<- unname(exp(rep(c(log(1), log(1), theta["log.rr.screen.m3"],log(1),log(1), log(1), theta["log.rr.screen.f3"]), each=2 ))) #screening rr by subpop
  rr.screen <- unname(exp(theta["log.rr.screen.ac"])) #screening rr in high activity group
  screen <- unname(t(apply(t(apply(screen_mat, 1, function(x) x*rr.screen.s)), 1, FUN=screen_fun, rr.screen=rr.screen))) #calculate true screening rates for each group by multiplying base rates by rr

  if (only_increasing_screening())
    screen <- apply(screen, 2, make_monotonically_increasing)

  screen[screen > 1] <- 1
  screen[screen < 0] <- 0

  # now we construct the column names for the screening matrix
  genders <-
    rep(c('male', 'female'), each=8)
  ethnicities <-
    rep(c('black', 'other', 'hispanic', 'msm'), each=2, times=2)
  ages <- rep(c('young', 'old'), times=8)
  pop_groups <- data.frame(sex=genders, demographic=ethnicities, age=ages)
  pop_groups <- dplyr::mutate(pop_groups[rep(seq_len(nrow(pop_groups)), each=2), ], risk=rep(c('low', 'high'), 16))
  pop_names <- apply(pop_groups, 1, paste, collapse=".")
  colnames(screen) <- pop_names
  return(screen)
}



#' @param K is carrying capacity
#' @param P0 is initial population
#' @param k is growth rate
#' @param t is a time vector
logistic_growth <- function(P0, K, k, t) {
  # K is the carrying capacity, as a percentage increase from P0.
  # This forces our estimates of screening to only go up.
  if (K == 0) return(rep(P0, length(t)))
  K <- (1-P0)*K + P0
  t <- t - mean(t)
  A = (K - P0)/P0
  logistic_curve <- K/(1+A*exp(-k*t))
  logistic_curve <- logistic_curve - min(logistic_curve)
  logistic_curve <- logistic_curve * (K-P0) / max(logistic_curve)
  logistic_curve <- logistic_curve + P0
  return(logistic_curve)
}


#' @param v a numeric vector to be made monotonically increasing.
make_monotonically_increasing <- function(v) {
  stopifnot(length(v) >= 2)
  stopifnot(is.numeric(v))

  for (i in 2:length(v))
    if (v[[i]] <= v[[i-1]]) v[[i]] <- v[[i-1]]

  return(v)
}

make_reporting_trend <- function(theta, ctrl.pt) {
  rep.bc<-ctrl.pt[[2]] #get internal control points for reporting bezier curve
  rep.bezier <- matrix(rep(c(ilogit(theta["logit.rep.symp.a"]), ilogit(theta["logit.rep.symp.d"]),rep.bc[1],rep.bc[2]),32), ncol=4, byrow=TRUE) #get 4 bezier points for calc reporting prob in each subgroup
  rep.trend <- apply(rep.bezier,1, bezier_fun) #baseline reporting probability over calibration period
  if (exists('gc_testing', envir=.GlobalEnv) && 'flat_screening' %in% gc_testing) {
    rep.bezier <- matrix(rep(rep(ilogit(theta["logit.rep.symp.a"]), 4),32), ncol=4, byrow=TRUE)
    rep.trend <- apply(rep.bezier,1, bezier_fun) #baseline reporting probability over calibration period
  }
  if (exists('gc_testing', envir=.GlobalEnv) && 'linear_screening' %in% gc_testing) {
    t_vec <-seq(0,1, length.out=(gc_env$cal.period+1))
    rep.trend <- matrix(rep(ilogit(theta['logit.rep.symp.a'])*(1-t_vec) + ilogit(theta['logit.rep.symp.d'])*t_vec, 32), ncol=32)
  }
  # if (logistic_screening_variation()) {
  #   t_vec <-seq(0,16, length.out=(gc_env$cal.period+1))
  #   rep.trend <- matrix(rep(c(ilogit(theta['logit.rep.symp.a']), ilogit(theta['logit.rep.symp.d']), exp(theta['log.rep.symp.growth'])), 32), ncol = 3, byrow = T)
  #   rep.trend <- apply(rep.trend, 1, function(x) logistic_growth(x[[1]], x[[2]], x[[3]], t_vec))
  # }
  # rep.trend <- apply(rep.trend, 2, make_monotonically_increasing)
  rep.trend[rep.trend>1]<-1  # make sure that reporting prob isn't >1
  return(rep.trend)
}
