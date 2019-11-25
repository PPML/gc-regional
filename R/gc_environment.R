#' The gcRegional Variable Environment
#'
#' @name gc_env
#'
#' @description This is the environment where variables and parameters necessary for the
#' gcRegional model will be stored. In prior versions of this model, many
#' variables were stored in the global scope -- any instances of those variables'
#' assignment and retrieval to the global scope have been replaced with assignment
#' and retrieval from gc_env.
#'
#' @section initial_files: The following gc_env objects are assigned in the
#'   initial_files function.
#'
#' \itemize{
#' \item s.dist.dat.m: Partnership distribution by subpop, Males, NSFG, not currently using for fitting
#' \item s.dist.dat.f: Partnership dist'n by subpop, Females, NSFG, not currently using for fitting
#' \item age.dist.dat: Proportion of partnerships with person of same age group, NSFG (national level estimates)
#' \item prev.extra.cat.dat: Prevalence, NHANES, not using for fitting but keep for comparison
#' \item diag.rr: RR of diagnosis rate for black and Hispanic, relative to overall rate for given age and sex
#' \item diag.age.sex.rate: Case rate data by age sex only, 2002-2016, sd calculated assuming +/-20%
#' \item p.symp.dat:
#' \item p.msm.dat: Proportion of male diagnosed cases in SSuN occuring in MSM, by age group
#' \item diag.subpop.rate: Rates by subpop for Baltimore, 2011-2016
#' \item priors: Prior distributions (parameters for probablility distributions for each input parameter)
#' \item theta: Vector of model parameter values
#' \item dat.s.dist: Proportion of within group mixing for each subpopulation, from NSFG data (not currently using for calibration)
#' \item var.s.dist.dat: Variance for mixing
#' \item s.dist.sd:
#' \item rr.diag.subpop: mean values of rr reported case by subpop
#' \item rr.diag.subpop.sd: sd for rr reported case by subpop (assuming +/- 20% and normal dist)
#' \item p.symp.ssun: Proportion of treated cases that are symptomatic
#' \item var.symp.ssun: Variance estimated from fitting data to beta distribution
#' \item p.msm.ssun: Proportion of males cases diagnosed in MSM, 2010-2016 SSuN
#' \item var.p.msm.ssun: Variance assuming beta distn and 95% CI=+/- 20% of mean
#' \item diag.rate: diagnosed case rates by age and sex only
#' \item diag.rate.sd: SD for diagnosed case data
#' \item prev.nhanes: GC prevalence 1999-2008 at national level, not currently using for model fitting
#' \item var.prev.nhanes:
#' \item sd.prior: starting value for standard deviation associated with each parameter, adapted during fitting
#' \item prior.param1: first parameter describing probablity distributions
#' \item prior.param2: second parameter describing probability distributions
#' }
#'
#' @section load_start: The following gc_env objects are assigned in the
#'   load_start function.
#'
#' \itemize{
#' \item site: The site location for modeling, either SF (San Francisco) or BA (Baltimore).
#' \item i: Number of subpopulations, 1= Black, 2=White, 3=Hispanic, 4=MS
#' \item j: Number of sexes, 1=male, 2=female
#' \item k:  Activity class, 1=low, 2=high
#' \item l: Age groups, 1= 15-24, 2=25-39
#' \item tstep: Weekly time step
#' \item cal.period:  Duration of calibration period (2002-2016 currently); Most data points are 2010 onwards, but allowing longer year period for time-varying parameters to prevent sudden changes
#' \item cal.start:  Time at which start time varying parameters are introduced
#' \item cal.end:  cal.start + cal.period
#' \item end.year:
#' \item start.year:
#' \item model.end:  For cpp code
#' \item age.cat:  age bandwidths, correspoonding to 15-24 and 25-39
#' \item out_all:
#' \item omega:  Supply and demand of sexual partnerships, 0.5=both sexes copmromise equally, 1=females determine number of partnerhips, 0=males determine number of partnerships
#' \item omega.t: Parameter describing compromise in partner number across subpops --> assuming that smaller pop size determines number of partnerships, rows= subpop of F, col=subpop of M. In SF,  hispanic pop is larger than black so black determines partnerships; in BA black pop larger than hisp and other.
#' \item n.total:  Total population size
#' \item p.msm:  Proportion of population MSM (x% of males) 4% estimate for state of Maryland; 18.5% for SF county in Grey et al. 2016
#' \item p.s.1:  Proportion of population in i=1 (black), Baltimore and SF 2016
#' \item p.s.3: Proportion of population in i=3 (hispanic), Baltimore and SF 2016
#' \item p.s.2: Proportion of population in i=2 (white), calculated as 1 - p.s.1 - p.s.3
#' \item p.sa.m.y:  Proportion of young M that are sexually active, by subpop (based on NSFG) - assuming mean for MSM
#' \item p.sa.m.o: Proportion of M old that are sexually active (based on NSFG)
#' \item p.sa.f.y: Proportion of young F that are sexually active (from NSFG, no sig diff across subpop)
#' \item p.sa.f.o:  Proportion of old F that are sexually active (from NSFG, no sig diff across subpop)
#' \item aging.rate.m: The rate at which sexually active males move from the young to age subpopulations.
#' \item aging.rate.f: The rate at which sexually active females move from the young to age subpopulations.
#' \item p.low: Proportion of the population in the low activity group
#' \item p.low.msm: Proportion of MSM in low activity group
#' \item p.sa.m: For each subpopulation, p.sa for M (y,o)
#' \item p.sa.f: For each subpopulation, p.sa for F (y,o)
#' \item p.sa: The proportion of each subpopulation (M/F, young/old, sexually-active/not-sexually-active)
#' \item ind: Subpopulation indexing matrix (k,l,j,i)
#' \item init.Y: Vector of initial number of infected, divide infected by symptomatic (Y) and asymptomatic (Z)
#' \item yinit: Initial distribution of the model population by sex, subpopulation, AC (M pop1 age1 L/H; M pop1 age 2 L/H; M pop2 age 1 L/H, M pop 2 age 2 L/H etc.)
#' \item s.index: Susceptible
#' \item y.index: Symptomatic infectious
#' \item z.index: Asymptomatic infectious
#' \item nsa.index: Not sexually active
#' \item diag.index: Diagnosed infections
#' \item symp.index: Diagnosed cases that are detected because symptomatic (rest are detected via screening)
#' \item inc.index: Incident infections
#' \item pop1: Subpop1
#' \item pop2: Subpop2
#' \item pop3: Subpop3
#' \item pop4: Subpop4
#' \item males:
#' \item females:
#' \item m1: M subpop1
#' \item m2: M subpop2
#' \item m3: M subpop3
#' \item m4: M subpop4
#' \item msw: M heterosexual
#' \item f1: F subpop1
#' \item f2: F subpop2
#' \item f3: F subpop3
#' \item f4: F subpop4
#' \item y.m: youngest age cat M
#' \item y.m.msw: yougest age cat MSW only
#' \item y.f: youngest age cat F
#' \item o.m: oldest age cat M
#' \item o.m.msw: oldest age cat MSW only
#' \item o.f: oldest age cat F
#' \item y.m.1: youngest age cat, subpop1
#' \item o.m.1: old age cat, subpop1
#' \item y.m.2: young age cat, subpop 2
#' \item o.m.2:  old age cat, subpop 2
#' \item y.m.3:  young age cat, subpop 3
#' \item o.m.3:  old age cat, subpop 3
#' \item y.m.4:  young age cat, subpop 4
#' \item o.m.4:  old age cat, subpop 4
#' \item y.f.1: youngest age cat, subpop1
#' \item o.f.1: old age cat, subpop1
#' \item y.f.2: young age cat, subpop 2
#' \item o.f.2:  old age cat, subpop 2
#' \item y.f.3:  young age cat, subpop 3
#' \item o.f.3:  old age cat, subpop 3
#' \item y.f.23:  young age cat, subpops 2&3
#' \item o.f.23:  old age cat, subpops 2&3
#' \item n.s.a:  pop sizes for given age, sex, subpop (sexually active only)
#' \item n.ns.a:  pop sizes for given age, sex, subpop (all)
#'
#'
#' @section pop_calc: The following gc_env objects are assigned in the
#'   pop_calc function.
#'
#' \itemize{
#' \item p.k: Distribution of population by low and high activity status
#' \item p.sa.array: Proportion of sexually active
#' \item n.dist: Population size
#' \item n.dist.sa: Population size for sexually active
#' \item p.dist: Population distribution within subpopulation and age group
#' \item p.s.dist: Relative sizes of different subpopulations - assume that mixing outside of subpopulation is proportionate to number of individuals in each subpopulation, excluding MSM
#' \item n.s.dist: Population sizes by sex (j) and subpopulation (i)
#' \item n.s.dist.sa: Population sizes by sex (j) and subpopulation (i) for sexually active population
#' \item n.s.pop: Population sizes, minus MSM
#' \item n.s.pop.sa: Population sizes, minus MSM, for sexually active population
#' \item n.i: (Population size by sex, subpopulation, AC, age) M, subpop=1:  age=AC=1, age=1 AC=2, age=2 AC=2; M subpop=2: age=AC=1,etc. Then F, subpop1...
#' \item n.sa: Size of sexually active population
#' \item n.nsa: Size of not sexually active population
#' }
#' @export
gc_env <- new.env(parent = globalenv())

#' Easier Named Assignment into the gc_env Environment
#'
#' The end result of this function is no different than
#' using `gc_env$var <- var` directly. Its utility is that it eliminates
#' the need to type the variables' names twice.
#' @export
gc_assign <- function(var) {
    invisible(assign(deparse(substitute(var)), var, envir=gc_env))
}

#' Easier Variable Retrieval from the gc_env Environment
#'
#' The functionality of gc_get is no different than performing
#' `varname <- gc_env$varname`. This function keeps the user
#' from needing to type the varname twice.
#' @export
gc_get <- function(varname, env = parent.env(environment())) {
  value <- get(varname, envir=gc_env)
  assign(varname, value, env=env)
}

