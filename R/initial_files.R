
#' Load Data from Initial Files
#'
#' This reads the data in the package's inst/extdata folder and
#' places that data into the gc_env environment.
#' @export
initial_files <- function() {

  #partnership dist'n by subpop, Males, NSFG, not currently using for fitting
  gc_env$s.dist.dat.m <-
    as.matrix(read.delim(
      system.file("extdata", "subpop_m_nsfg.txt", package = 'gcRegional'),
      row.names = 1
    ))

  #partnership dist'n by subpop, Females, NSFG, not currently using for fitting
  # s.dist.dat.f <- as.matrix(read.delim("R_inputs/subpop_f_nsfg.txt", row.names=1))
  gc_env$s.dist.dat.f <-
    as.matrix(read.delim(
      system.file("extdata", "subpop_f_nsfg.txt", package = 'gcRegional'),
      row.names = 1
    ))

  #proportion of partnerships with person of same age group, NSFG (national level estimates)
  # age.dist.dat <- read.delim("R_inputs/age_nsfg.txt")
  gc_env$age.dist.dat <-
    read.delim(system.file("extdata", "age_nsfg.txt", package = 'gcRegional'))
  popsizes <- estBetaParamsVar(gc_env$age.dist.dat$p.same.age, gc_env$age.dist.dat$var)
  gc_env$age.dist.dat$N <- popsizes[[1]] + popsizes[[2]]

  #prevalence, NHANES, not using for fitting but keep for comparison
  # prev.extra.cat.dat <- read.delim("R_inputs/prev_extra_cat.txt")
  gc_env$prev.extra.cat.dat <-
    read.delim(system.file("extdata", "prev_extra_cat.txt", package = 'gcRegional'))

  # NHANES prevalence data
  gc_env$nhanes.updated.dat <-
    readRDS(system.file("data", "nhanes_prevalence_99-08.rds", package = 'gcRegional'))

  #RR of diagnosis rate for black and Hispanic, relative to overall rate for given age and sex
  # diag.rr <- read.delim(paste("R_inputs/subpop_rr_", gc_env$site, ".txt", sep=""))
  gc_env$diag.rr <-
    read.delim(system.file(
      "extdata",
      paste("subpop_rr_", gc_env$site, ".txt", sep = ""),
      package = 'gcRegional'
    ))

  #case rate data by age sex only, 2002-2016, sd calculated assuming +/-20%
  # diag.age.sex.rate <- as.matrix(read.delim(paste("R_inputs/diag_age_sex_rate_", gc_env$site,".txt", sep=""), row.names=1))
  gc_env$diag.age.sex.rate <-
    as.matrix(read.delim(
      system.file(
        "extdata",
        paste("diag_age_sex_rate_", gc_env$site, ".txt", sep = ""),
        package = 'gcRegional'
      ),
      row.names = 1
    ))


  # p.symp.dat <- read.delim(paste("R_inputs/ssun_trt_sympt_", gc_env$site, ".txt", sep=""))
  gc_env$p.symp.dat <-
    read.delim(system.file(
      "extdata",
      paste("ssun_trt_sympt_", gc_env$site, ".txt", sep = ""),
      package = 'gcRegional'
    ))

  #proportion of male diagnosed cases in SSuN occuring in MSM, by age group
  # p.msm.dat <- read.delim(paste("R_inputs/ssun_p_msm_", gc_env$site, ".txt", sep=""))
  gc_env$p.msm.dat <-
    read.delim(system.file(
      "extdata",
      paste("ssun_p_msm_", gc_env$site, ".txt", sep = ""),
      package = 'gcRegional'
    ))

  #Rates by subpop for Baltimore, 2011-2016
  # diag.subpop.rate <-  as.data.frame(read.delim(paste("R_inputs/diag_subpop_rate_", gc_env$site, ".txt", sep=""), row.names=1))
  gc_env$diag.subpop.rate <-
    as.data.frame(read.delim(
      system.file(
        "extdata",
        paste("diag_subpop_rate_", gc_env$site, ".txt", sep = ""),
        package = 'gcRegional'
      ),
      row.names = 1
    ))

  gc_env$diag.subpop.rate.denom <-
    readRDS(system.file("extdata", paste0(gc_env$site, "_subpop_rate_denom.rds"), package='gcRegional'))

  # load starting value of parmameters
  # load(paste("R_inputs/theta_", gc_env$site, ".rda", sep=""))
  load(system.file(
    "extdata",
    paste("theta_", gc_env$site, ".rda", sep = ""),
    package = 'gcRegional'
  ))
  gc_env$theta <- theta

  # theta <- readRDS(system.file("extdata", "theta_with_national_means.rds", package='gcRegional'))

  # get national natural history parameters
  gc_env$national_natural_history <- readRDS(system.file("data", "national_natural_history_parameters.rds", package = 'gcRegional'))

  ### prep data used for calibration

  #proportion of within group mixing for each subpopulation, from NSFG data (not currently using for calibration)
  gc_env$dat.s.dist <-
    c(diag(gc_env$s.dist.dat.m), diag(gc_env$s.dist.dat.f))

  #variance for mixing
  gc_env$var.s.dist.dat <-
    c(gc_env$s.dist.dat.m[, "var"], gc_env$s.dist.dat.f[, "var"])
  gc_env$s.dist.sd <-
    as.data.frame(
      cbind(
        gc_env$dat.s.dist,
        rbind(gc_env$s.dist.dat.m[, 4:5], gc_env$s.dist.dat.f[, 4:5])
      ),
      row.names = c(
        "black M",
        "other M",
        "hispanic M",
        "black F",
        "other F",
        "hispanic F"
      )
    )

  #mean values of rr reported case by subpop
  gc_env$rr.diag.subpop <- as.numeric(gc_env$diag.rr[, "rr_diag"])

  #sd for rr reported case by subpop (assuming +/- 20% and normal dist)
  gc_env$rr.diag.subpop.sd <- as.numeric(gc_env$diag.rr[, "sd"])

  #proportion of treated cases that are symptomatic
  gc_env$p.symp.ssun <- gc_env$p.symp.dat$p.symp

  #variance estimated from fitting data to beta distribution
  gc_env$var.symp.ssun <- gc_env$p.symp.dat$p.symp.var

  #proportion of males cases diagnosed in MSM, 2010-2016 SSuN
  gc_env$p.msm.ssun <- gc_env$p.msm.dat$p_msm

  #variance assuming beta distn and 95% CI=+/- 20% of mean
  gc_env$var.p.msm.ssun <- gc_env$p.msm.dat$var

  #diagnosed case rates by age and sex only
  gc_env$diag.rate <-
    as.vector(gc_env$diag.age.sex.rate[, c("m_y", "m_o", "f_y", "f_o")])

  #SD for diagnosed case data for 2002-2016
  gc_env$diag.rate.sd <-
    as.vector(gc_env$diag.age.sex.rate[, c("m_y_sd", "m_o_sd", "f_y_sd", "f_o_sd")])

  #GC prevalence 1999-2008 at national level, not currently using for model fitting
  gc_env$prev.nhanes <- as.numeric(gc_env$prev.extra.cat.dat$prev)
  gc_env$var.prev.nhanes <- as.numeric(gc_env$prev.extra.cat.dat$var)

  gc_env$p.symp.popsizes <-
    readRDS(
      system.file(
        "extdata",
        paste0(gc_env$site, "_p_symp_popsizes.rds"),
        package = 'gcRegional'))


  ### read in model priors (used by prior.function.R)
  load_priors()
}


load_priors <- function() {

  # prior distributions (parameters for probablility distributions for each input parameter)
  # priors <- read.delim("R_inputs/priors.txt")
  gc_env$priors <-
    read.delim(system.file("extdata", "priors.txt", package = 'gcRegional'))

  #starting value for standard deviation associated with each parameter, adapted during fitting
  gc_env$sd.prior <- gc_env$priors$sd.transf.1

  #first parameter describing probablity distribution
  gc_env$prior.param1 <-  gc_env$priors$param1
  prior.param1 <- gc_env$priors$param1
  names(gc_env$prior.param1) <- gc_env$priors$parameter

  #second parameter describing probability distributions
  gc_env$prior.param2 <-  gc_env$priors$param2
  names(gc_env$prior.param2) <- gc_env$priors$parameter
}
