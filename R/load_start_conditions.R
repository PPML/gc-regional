#' Load Starting Conditions
#'
#' When running load_start() the user will be prompted to specify the site location.
#' After doing so, initial_files() will run in order to load the initial parameters
#' and variables that correspond to the specified site. The site and variables
#' loaded will be stored in the gc_env environment.
#' @export
load_start <- function(site=NA) {

    # A local function to help assign the gc_env$site variable.
    assign_site <- function(site) {
        if (site %in% c("SF", "BA")) gc_env$site <- site
        else gc_env$site <- readline("Enter the site (SF for San Francisco or BA for Baltimore): ")
    }

    # Ask the user to update gc_env$site until they answer with either SF or BA.
    assign_site(site)
    while(length(gc_env$site) == 0 || !gc_env$site %in% c("SF", "BA")) {
        if (length(gc_env$site) != 0) cat("The site must be one of SF or BA.\n")
        assign_site(site)
    }
    site <- gc_env$site

    initial_files() # load required data files
    load_SF16_data()
    load_BA17_data()
    load_BA16_data()

  ### set up different subpopulations, sexes, and activity classes ###
  i<-4 #number of subpopulations, 1= Black, 2=White, 3=Hispanic, 4=MSM
  j<-2 # number of sexes, 1=male, 2=female
  k<-2 #activity class, 1=low, 2=high
  l<-2 #age groups, 1= 15-24; 2=25-39
  index <- i*j*k*l
  gc_assign(i)
  gc_assign(j)
  gc_assign(k)
  gc_assign(l)
  gc_assign(index)

  #### model time steps and calibration period ####
  tstep <- 1/52 # weekly time step

  #duration of calibration period (2002-2016 currently); most data points are
  #2010 onwards, but allowing longer year period for time-varying parameters to
  #prevent sudden changes
  cal.period <- ifelse(site == 'BA', 16, 15) #duration of calibration period (2002-2016 currently); most data points are 2010 onwards, but allowing longer year period for time-varying parameters to prevent sudden changes
  cal.start <- 60  #time at which start introducing time varying parameters
  cal.end <- cal.start + cal.period
  end.year <- ifelse(site == 'BA', 2017, 2016)
  start.year <- 2002
  intervention_period <- 4
  model.end <- cal.start + cal.period + intervention_period + 1 #for cpp code
  cal_and_intervention_rows <- (2010-1942+1):(model.end)

  gc_assign(tstep)
  gc_assign(cal.period)
  gc_assign(cal.start)
  gc_assign(cal.end)
  gc_assign(end.year)
  gc_assign(start.year)
  gc_assign(model.end)
  gc_assign(intervention_period)
  gc_assign(cal_and_intervention_rows)

  age.cat<-c(10, 15) # age band widths, corresponding to 15-24 yo and 25-39 yo
  out_all <- NULL

  #supply and demand of sexual partnerships, 0.5=both sexes copmromise equally,
  #1=females determine number of partnerhips, 0=males determine number of
  #partnerships
  omega <- 0.5

  #parameter describing compromise in partner number across subpops --> assuming
  #that smaller pop size determines number of partnerships, rows= subpop of F,
  #col=subpop of M
  omega.t <-
    if (site == "BA")
      matrix(c(0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0),
             nrow = 4,
             byrow = T)
  else
    matrix(c(0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0),
           nrow = 4,
           byrow = T)
  #in SF, hispanic pop is larger than black, so black determines partnerships;
  #in BA, black pop larger than hisp and other
  gc_assign(age.cat)
  gc_assign(out_all)
  gc_assign(omega)
  gc_assign(omega.t)


  #### population characteristics for Baltimore and San Francisco ###

  # total population size, Baltimore or SF 2016
  n.total <- if(site=="BA") 242185 else 359020

  # proportion of population MSM (x% of males) 4% estimate for state of
  # Maryland; 18.5% estimate for SF county  in Grey et al. 2016
  p.msm <- if(site=="BA") 0.04 else 0.185
  if (site=="SF" && lower_proportion_msm()) p.msm <- 0.09

  #proportion of population in i=1 (black), Baltimore and SF 2016
  p.s.1 <-if(site=="BA") 0.57 else 0.045

  #proportion of population in i=3 (hispanic), Baltimore and SF 2016
  p.s.3 <-if(site=="BA") 0.06 else 0.162
  p.s.2<-1-p.s.1-p.s.3

  #proportion of young M that are sexually active, by subpop (based on NSFG) -
  # assuming mean for MSM
  p.sa.m.y <- c(0.78, 0.64, 0.68, 0.67)

  #proportion of M old that are sexually active (based on NSFG)
  p.sa.m.o <-rep(0.96,4)

  #proportion of young  F that are sexually active (from NSFG, no sig diff across subpop)
  p.sa.f.y <-rep(0.66,4)

  #proportion of old F that are sexually active (from NSFG, no sig diff across subpop)
  p.sa.f.o <- rep(0.98,4)
  aging.rate.m <- (p.sa.m.o - p.sa.m.y) / (1-p.sa.m.y)
  aging.rate.f <- (p.sa.f.o - p.sa.f.y) / (1-p.sa.f.y)

  #proportion of population in low activity group
  p.low.cond <- exists("p.low", envir=.GlobalEnv)
  p.low <- ifelse(p.low.cond, get("p.low", envir = .GlobalEnv), 0.90)
  p.low.msm <- ifelse(p.low.cond, get("p.low", envir = .GlobalEnv), 0.90) # proportion of MSM in low activity group

  #for each subpop, p.sa for M (y,o)
  p.sa.m<-array(c(p.sa.m.y, p.sa.m.o),dim=c(4,2))

  #for each subpop, p.sa for M (y,o)
  p.sa.f<-array(c(p.sa.f.y, p.sa.f.o),dim=c(4,2))
  p.sa<-c(as.vector(rep(aperm(p.sa.m, perm=c(2,1)), each=2)), as.vector(rep(aperm(p.sa.f, perm=c(2,1)), each=2)))

  gc_assign(n.total)
  gc_assign(p.msm)
  gc_assign(p.s.1)
  gc_assign(p.s.3)
  gc_assign(p.s.2)
  gc_assign(p.sa.m.y)
  gc_assign(p.sa.m.o)
  gc_assign(p.sa.f.y)
  gc_assign(p.sa.f.o)
  gc_assign(aging.rate.m)
  gc_assign(aging.rate.f)
  gc_assign(p.low )
  gc_assign(p.low.msm)
  gc_assign(p.sa.m)
  gc_assign(p.sa.f)
  gc_assign(p.sa)


  #for each subpopulation (i), need to define pop. size (M+F) and distribution of activity classes and ages
  pop_calc(n.total)
  aging_fun(age.cat)

  ### initial model conditions ###
  ind<- as.matrix(expand.grid(1:k,1:l,1:j,1:i)) #indexing matrix (k,l,j,i)

  #vector of initial number of infecteds,  divide infecteds by symptomatic (Y) and asymptomatic (Z)
  init.Y <-c(rep(c(20,20,20,20),4),rep(c(20,20,20,20),3),c(0,0,0,0))
  gc_assign(ind)
  gc_assign(init.Y)

  # initial distribution of the model population by sex, subpop, AC (M pop1
  # age1 L/H; M pop1 age 2 L/H; M pop2 age 1 L/H, M pop 2 age 2 L/H, etc.)
  n.sa <- gc_env$n.sa
  n.nsa <- gc_env$n.nsa
  n.i <- gc_env$n.i

  yinit <- c(
    S=n.sa-init.Y, #susceptible
    Y=init.Y*0.5, #infectious, symptomatic
    Z=init.Y*0.5, #infectious, asymptomatic
    NSA=n.nsa,  #not sexually active
    D=rep(0,index), #diagnosed cases (cumulative)
    R=rep(0,index), #diagnosed symptomatic cases (cumulative)
    CUMINC = rep(0,index), # cumulative incidence
    SCR = rep(0,index) # cumulative number of tests
  )
  gc_assign(yinit)

  #indexes for pulling out different states from the model output matrix
  s.index <- 1:index  #susceptible
  y.index <- index+1:index #symptomatic infectious
  z.index <- index*2+1:index #asymptomatic infectious
  nsa.index <- index*3+1:index #not sexually active
  diag.index <- index*4+1:index #diagnosed infections
  symp.index <- index*5+1:index #diagnosed cases that are detected because symptomatic (rest are detected via screening)
  inc.index <- index*6+1:index #incident infections
  scr.index <- index*7+1:index

  gc_assign(s.index)
  gc_assign(y.index)
  gc_assign(z.index)
  gc_assign(nsa.index)
  gc_assign(diag.index)
  gc_assign(symp.index)
  gc_assign(inc.index)
  gc_assign(scr.index)

  pop1 <- c(1:4,17:20) #subpop1
  pop2 <- c(5:8,21:24) # subpop2
  pop3 <- c(9:12,25:28) #subpop3
  pop4 <-c(13:16,29:32 ) #subpop4
  males <-1:16 # males
  females <-17:32 #females
  m1<-1:4 #M subpop1
  m2<-5:8 #M subpop2
  m3<-9:12 #M subpop3
  m4<-13:16 #M subpop4
  msw <-1:12 #M heterosexual
  f1<-17:20 #F subpop1
  f2<-21:24 #F subpop2
  f3<-25:28 #F subpop3
  f4<-29:32 #F subpop4
  y.m<-c(1:2,5:6,9:10,13:14) #youngest age cat M
  y.m.msw<- c(1:2, 5:6, 9:10) #yougest age cat MSW only
  y.f<-c(17:18,21:22,25:26) #youngest age cat F
  o.m<-c(3:4,7:8, 11:12,15:16) #oldest age cat M
  o.m.msw <- c(3:4, 7:8, 11:12) #oldest age cat MSW only
  o.f<-c(19:20,23:24,27:28) #oldest age cat F
  y.m.1 <- 1:2 #youngest age cat, subpop1
  o.m.1 <- 3:4 #old age cat, subpop1
  y.m.2 <- 5:6 #young age cat, subpop 2
  o.m.2 <- 7:8 # old age cat, subpop 2
  y.m.3 <- 9:10 # young age cat, subpop 3
  o.m.3 <- 11:12 # old age cat, subpop 3
  y.m.4 <- 13:14 # young age cat, subpop 4
  o.m.4 <- 15:16 # old age cat, subpop 4
  y.f.1 <- 17:18 #youngest age cat, subpop1
  o.f.1 <- 19:20 #old age cat, subpop1
  y.f.2 <- 21:22 #young age cat, subpop 2
  o.f.2 <- 23:24 # old age cat, subpop 2
  y.f.3 <- 25:26 # young age cat, subpop 3
  o.f.3 <- 27:28 # old age cat, subpop 3
  y.f.23 <- c(21:22,25:26) # young age cat, subpops 2&3
  o.f.23 <- c(23:24,27:28) # old age cat, subpops 2&3

  # pop sizes for given age, sex, subpop (sexually active only)
  n.s.a<-(sapply(1:((length(diag.index)-4)/2),function(x){sum(n.sa[1:2+(x-1)*2])}))

  # pop sizes for given age, sex, subpop (all)
  n.ns.a<-(sapply(1:((length(diag.index)-4)/2),function(x){sum(n.i[1:2+(x-1)*2])}))

  gc_assign(pop1)
  gc_assign(pop2)
  gc_assign(pop3)
  gc_assign(pop4)
  gc_assign(males)
  gc_assign(females)
  gc_assign(m1)
  gc_assign(m2)
  gc_assign(m3)
  gc_assign(m4)
  gc_assign(msw)
  gc_assign(f1)
  gc_assign(f2)
  gc_assign(f3)
  gc_assign(f4)
  gc_assign(y.m)
  gc_assign(y.m.msw)
  gc_assign(y.f)
  gc_assign(o.m)
  gc_assign(o.m.msw)
  gc_assign(o.f)
  gc_assign(y.m.1)
  gc_assign(o.m.1)
  gc_assign(y.m.2)
  gc_assign(o.m.2)
  gc_assign(y.m.3)
  gc_assign(o.m.3)
  gc_assign(y.m.4)
  gc_assign(o.m.4)
  gc_assign(y.f.1)
  gc_assign(o.f.1)
  gc_assign(y.f.2)
  gc_assign(o.f.2)
  gc_assign(y.f.3)
  gc_assign(o.f.3)
  gc_assign(y.f.23)
  gc_assign(o.f.23)
  gc_assign(n.s.a)
  gc_assign(n.ns.a)


  intervention_names <-
    c(
      base_case = "Base Case",
      no_screening = "No Screening",
      universal_annual = "Universal Annual",
      universal_semi_annual = "Universal 2x Annual",
      high_activity_young_annual = "High Activity Young Annual Screening",
      high_activity_young_semi_annual = "High Activity Young 2x Annual Screening",
      high_activity_young_quarter_annual = "High Activity Young 4x Annual Screening",
      young_annual = "Young Annual Screening",
      young_semi_annual = "Young 2x Annual Screening",
      young_quarter_annual = "Young 4x Annual Screening",
      high_activity_MSM_annual = "High Activity MSM Annual Screening",
      high_activity_MSM_semi_annual = "High Activity MSM 2x Annual Screening",
      high_activity_MSM_quarter_annual = "High Activity MSM 4x Annual Screening",
      MSM_annual = "MSM Annual Screening",
      MSM_semi_annual = "MSM 2x Annual Screening",
      MSM_quarter_annual = "MSM 4x Annual Screening",
      high_activity_annual = "High Activity Annual Screening",
      high_activity_semi_annual = "High Activity 2x Annual Screening",
      high_activity_quarter_annual = "High Activity 4x Annual Screening",
      female_young_annual = "Female Young Annual Screening",
      female_young_semi_annual = "Female Young 2x Annual Screening",
      female_young_quarter_annual = "Female Young 4x Annual Screening",
      remove_ltfu_10pct = "Remove 10% LTFU",
      remove_ltfu_20pct = "Remove 20% LTFU",
      remove_f_msm_10pct_and_msw_20pct_ltfu = "Remove 10% LTFU for Women and MSM, 20% LTFU for MSW",
      `1year_high_activity_semi_annual` = 'One Year Van Intervention'
    )

  gc_assign(intervention_names)

  sd.prior <- gc_env$sd.prior

  theta <- gc_env$theta

  if (logistic_screening_variation()) {
    screening_names <- c(
      "m1.1", "m2.1", "m1", "m2", "m1", "m2", 'msm1', 'msm2', "f1.1", "f2.1", "f1", "f2", "f1", "f2")
    for (scr_name in unique(screening_names)) {
      theta[paste0('log.', scr_name, '.growth')] <- log(0.2)
      sd.prior[length(sd.prior)+1] <- exp(sqrt(2/9))
    }
    gc_env$screening_names <- screening_names
    theta['log.rep.symp.growth'] <- log(0.2)
    sd.prior[length(sd.prior)+1] <- exp(log(sqrt(2/9)))

    gc_env$theta <- theta
    gc_env$sd.prior <- sd.prior
    adjust_prior_for_logistic_screening()
  }

  if (check_variation('population_risk_mixing')) {
    adjust_priors_for_population_risk_mixing()
    gc_env$theta['logit.epsilon'] <- logit(rbeta(1, 1.1, 1.1))
  }

  construct_key_populations()

  gc_env$quantile_funs <- list(
    beta = qbeta,
    gamma = qgamma,
    normal = qnorm
  )

  gc_env$transforms <- list(
    log = log,
    logit = logit
  )

  gc_env$random_funs <- list(
    beta = rbeta,
    gamma = rgamma,
    normal = rnorm
  )

  gc_env$distribution_density_lookup <- list(
    beta = dbeta,
    gamma = dgamma,
    normal = dnorm
  )

  gc_env$inverse_transforms <- list(
    logit = ilogit,
    log = exp
  )

  if (check_variation('national_posterior_medians')) {
    gc_env$national_posterior_medians <- apply(gc_env$national_natural_history, 2, median)
    national_params <-
      sapply(names(gc_env$national_posterior_medians), function(x) {
        x <-
          strsplit(x, "\\.")
        paste0(x[[1]][2:length(x[[1]])], collapse = '.')
      })
    gc_env$priors <- gc_env$priors[! gc_env$priors$parameter %in% national_params, ]
    gc_env$theta <- remove_national_parameters(gc_env$theta)
  }

  if (check_variation('no_differential_reporting')) {
    # remove differential reporting probabilities for symptomatic cases
    for_removal <- grepl("risk.rep", gc_env$priors$parameter)
    gc_env$priors <- gc_env$priors[! for_removal, ]
    gc_env$theta <- gc_env$theta[! for_removal]
  }


  if (check_variation('widen_race_assortativity')) {
    gc_env$priors[gc_env$priors$parameter %in% c('theta.1', 'theta.2', 'theta.3', 'theta.5', 'theta.6', 'theta.7'), 'param1'] <- 5.833333
    gc_env$priors[gc_env$priors$parameter %in% c('theta.1', 'theta.2', 'theta.3', 'theta.5', 'theta.6', 'theta.7'), 'param2'] <- 2.5
    gc_env$priors[gc_env$priors$parameter == 'theta.4', 'param1'] <- 9.866914
    gc_env$priors[gc_env$priors$parameter == 'theta.4', 'param2'] <- 1.219506
  }

  if (check_variation('widen_beziers')) {
    gc_env$priors[grepl("screen", gc_env$priors$parameter), c('param1', 'param2')] <- 1.1
  }

  if (check_variation('separate_increased_risk')) {
    gc_env$theta['logit.behav.hetero'] <- logit(0)
    gc_env$priors <- rbind.data.frame(
      gc_env$priors,
      data.frame(
        parameter = 'behav.hetero',
        distribution = 'beta',
        param1 = 1,
        param2 = 15,
        mean = .5,
        sd = 0,
        sd.prior = .59,
        transformation = 'logit',
        mean.transf = 0,
        ucl.transf = 0,
        sd.transf.1 = 0,
        sd.transf = 0,
        sd1 = 0
      )
    )

    #first parameter describing probablity distribution
    gc_env$prior.param1 <-  gc_env$priors$param1
    names(gc_env$prior.param1) <- gc_env$priors$parameter

    #second parameter describing probability distributions
    gc_env$prior.param2 <-  gc_env$priors$param2
    names(gc_env$prior.param2) <- gc_env$priors$parameter

  }

}

