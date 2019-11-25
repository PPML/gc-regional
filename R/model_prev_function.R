######################################################
### function to calculate outputs used for fitting ###
######################################################

model_prev <- function(sol) {  #sol is the output from the transmission model, a matrix of model compartment sizes at each time step

  # Retrieve Environment Variables
  s.index <- gc_env$s.index
  y.index <- gc_env$y.index
  z.index <- gc_env$z.index
  nsa.index <- gc_env$nsa.index
  cal.start <- gc_env$cal.start
  cal.end <- gc_env$cal.end
  msw <- gc_env$msw
  females <- gc_env$females
  males <- gc_env$males
  y.m <- gc_env$y.m
  y.f <- gc_env$y.f
  y.m.msw <- gc_env$y.m.msw
  y.m.4 <- gc_env$y.m.4
  o.f <- gc_env$o.f
  o.m <- gc_env$o.m
  pop1 <- gc_env$pop1
  pop2 <- gc_env$pop2
  pop3 <- gc_env$pop3
  pop4 <- gc_env$pop4
  y.f.1 <- gc_env$y.f.1
  y.f.2 <- gc_env$y.f.2
  y.f.3 <- gc_env$y.f.3
  y.f.23 <- gc_env$y.f.23
  o.f.1 <- gc_env$o.f.1
  o.f.2 <- gc_env$o.f.2
  o.f.3 <- gc_env$o.f.3
  o.f.23 <- gc_env$o.f.23
  y.m.1 <- gc_env$y.m.1
  y.m.2 <- gc_env$y.m.2
  y.m.3 <- gc_env$y.m.3
  o.m.1 <- gc_env$o.m.1
  o.m.2 <- gc_env$o.m.2
  o.m.3 <- gc_env$o.m.3
  o.m.4 <- gc_env$o.m.4
  inc.index <- gc_env$inc.index
  diag.age.sex.rate <- gc_env$diag.age.sex.rate
  diag.subpop.rate <- gc_env$diag.subpop.rate
  diag.index <- gc_env$diag.index
  diag.rr <- gc_env$diag.rr
  symp.index <- gc_env$symp.index
  o.m.msw <- gc_env$o.m.msw
  p.msm.dat <- gc_env$p.msm.dat
  age.dist.dat <- gc_env$age.dist.dat
  var.symp.ssun <- gc_env$var.symp.ssun
  diag.rate.sd <- gc_env$diag.rate.sd
  var.p.msm.ssun <- gc_env$var.p.msm.ssun
  m1 <- gc_env$m1
  m2 <- gc_env$m2
  m3 <- gc_env$m3
  m4 <- gc_env$m4
  n.i <- gc_env$n.i
  f1 <- gc_env$f1
  f2 <- gc_env$f2
  f3 <- gc_env$f3
  f4 <- gc_env$f4
  p.s.1 <- gc_env$p.s.1
  p.s.2 <- gc_env$p.s.2
  p.s.3 <- gc_env$p.s.3
  p.s.4 <- gc_env$p.s.4
  cal.period <- gc_env$cal.period

  start_year <- gc_env$cal.start + 8
  end_year <- gc_env$cal.end


  #calculate total population size, sexually active (sa) and total
  pop.size.sa <- sol[,1+s.index] + sol[,1+y.index] + sol[,1+z.index]
  pop.size.all <- sol[,1+s.index] + sol[,1+y.index] + sol[,1+z.index] + sol[,1+nsa.index]
  #calculate total prevalent infections, and by different stratifications (to align with NHANES estimates if comparing with national prevalence estimates)
  #adjusted values are redistributing MSM cases among the race/ethnic subpopulations
  inf.total<-sol[,1+y.index]+sol[,1+z.index]
  inf.subpop.age <-colSums(matrix(tail(inf.total,1), nrow=2))
  pop.subpop.age <- colSums(matrix(tail(pop.size.all,1), nrow=2))
  pop.sub.adj <- rep(pop.subpop.age[9:14],2) #pop size by sex, subpop, redistributing MSM for calculating reported rates
  pop.msm <- pop.subpop.age[7:8]
  prev.total <- cbind(rowSums(inf.total[,m1])/sum(n.i[m1]), rowSums(inf.total[,m2])/sum(n.i[m2]),rowSums(inf.total[,m3])/sum(n.i[m3]),rowSums(inf.total[,m4])/sum(n.i[m4]),rowSums(inf.total[,f1])/sum(n.i[f1]),rowSums(inf.total[,f2])/sum(n.i[f2]),rowSums(inf.total[,f3])/sum(n.i[f3]))
  prev.age.subpop.m.adj.msm <- cbind((rowSums(inf.total[,m1[1:2]])+rowSums(inf.total[,m4[1:2]])*p.s.1)/(rowSums(pop.size.all[,m1[1:2]])+rowSums(pop.size.all[,m4[1:2]])*p.s.1),
                                          (rowSums(inf.total[,m1[1:4]])+rowSums(inf.total[,m4[1:4]])*p.s.1)/(rowSums(pop.size.all[,m1[1:4]])+rowSums(pop.size.all[,m4[1:4]])*p.s.1),
                                     (rowSums(inf.total[,m2[1:2]])+rowSums(inf.total[,m3[1:2]])+rowSums(inf.total[,m4[1:2]])*(p.s.2+p.s.3))/(rowSums(pop.size.all[,m2[1:2]])+rowSums(pop.size.all[,m3[1:2]])+rowSums(pop.size.all[,m4[1:2]])*(p.s.2+p.s.3)),
                                     (rowSums(inf.total[,m2[1:4]])+rowSums(inf.total[,m3[1:4]])+rowSums(inf.total[,m4[1:4]])*(p.s.2+p.s.3))/(rowSums(pop.size.all[,m2[1:4]])+rowSums(pop.size.all[,m3[1:4]])+rowSums(pop.size.all[,m4[3:4]])*(p.s.2+p.s.3))
                                     )
  prev.age.subpop.m.adj.msm<-colMeans(prev.age.subpop.m.adj.msm[(start_year):(end_year),])
  prev <- mean(rowSums(inf.total[(start_year):(end_year),])/rowSums(pop.size.all[(start_year):(end_year),])) #overall population prevalence
  prev.adj <- mean((rowSums(inf.total[(start_year):(end_year),msw]) + rowSums(inf.total[(start_year):(end_year),m4]) + rowSums(inf.total[(start_year):(end_year),females]))/rowSums(pop.size.all[(start_year):(end_year),]))
  prev.m <- mean(rowSums(inf.total[(start_year):(end_year),males])/rowSums(pop.size.all[(start_year):(end_year),males]))
  prev.m.adj <- mean((rowSums(inf.total[(start_year):(end_year),msw]) + rowSums(inf.total[(start_year):(end_year),m4]))/rowSums(pop.size.all[(start_year):(end_year),males]) )
  prev.f<- mean(rowSums(inf.total[(start_year):(end_year),females])/rowSums(pop.size.all[(start_year):(end_year),females]))
  prev.y <- mean((rowSums(inf.total[(start_year):(end_year),y.m]) + rowSums(inf.total[(start_year):(end_year),y.f]))/(rowSums(pop.size.all[(start_year):(end_year),y.m]) + rowSums(pop.size.all[(start_year):(end_year),y.f]))) #overall prevalence in age cat 1
  prev.y.adj <- mean((rowSums(inf.total[(start_year):(end_year),y.m.msw])+ rowSums(inf.total[(start_year):(end_year),y.m.4]) + rowSums(inf.total[(start_year):(end_year),y.f]))/(rowSums(pop.size.all[(start_year):(end_year),y.m]) + rowSums(pop.size.all[(start_year):(end_year),y.f]))) #overall prevalence in age cat 1
  prev.y.m <- mean(rowSums(inf.total[(start_year):(end_year),y.m])/rowSums(pop.size.all[(start_year):(end_year),y.m]))#prevalence in M, age cat 1
  prev.y.msw <- mean(rowSums(inf.total[(start_year):(end_year),y.m.msw])/rowSums(pop.size.all[(start_year):(end_year),y.m.msw])) #prevalence in MSW, age cat 1
  prev.y.m.adj <-mean((rowSums(inf.total[(start_year):(end_year),y.m.msw])+ rowSums(inf.total[(start_year):(end_year),y.m.4]))/rowSums(pop.size.all[(start_year):(end_year),y.m]))
  prev.y.f <- mean(rowSums(inf.total[(start_year):(end_year),y.f])/rowSums(pop.size.all[(start_year):(end_year),y.f])) #prevalence in F, age cat 1
  prev.o.m <- mean(rowSums(inf.total[(start_year):(end_year),o.m])/rowSums(pop.size.all[(start_year):(end_year),o.m])) #prevalence in M, age cat 2
  prev.o.f <- mean(rowSums(inf.total[(start_year):(end_year),o.f])/rowSums(pop.size.all[(start_year):(end_year),o.f])) #prevalence in F, age cat 2
  prev.1.adj <- mean((rowSums(inf.total[(start_year):(end_year),pop1])+rowSums(inf.total[(start_year):(end_year),pop4])*p.s.1)/(rowSums(pop.size.all[(start_year):(end_year),pop1])+ rowSums(pop.size.all[(start_year):(end_year),pop4])*p.s.1 ))#overall prevalence in subpop 1
  prev.1.no.msm <- mean(rowSums(inf.total[(start_year):(end_year),pop1])/rowSums(pop.size.all[(start_year):(end_year),pop1]))#overall prevalence in subpop 1
  prev.23.adj <- mean((rowSums(inf.total[(start_year):(end_year),pop2]) + rowSums(inf.total[(start_year):(end_year),pop3]) + rowSums(inf.total[(start_year):(end_year),pop4])*(1-p.s.1) ) / (rowSums(pop.size.all[(start_year):(end_year),pop2])+ rowSums(pop.size.all[(start_year):(end_year),pop3]) + rowSums(pop.size.all[(start_year):(end_year),pop4])*(1-p.s.1) )) #overall prevalence in subpops 2+3
  prev.23.no.msm <- mean((rowSums(inf.total[(start_year):(end_year),pop2]) + rowSums(inf.total[(start_year):(end_year),pop3])) / (rowSums(pop.size.all[(start_year):(end_year),pop2])+ rowSums(pop.size.all[(start_year):(end_year),pop3]) )) #overall prevalence in subpops 2+3
  prev.1.f <- mean(rowSums(inf.total[(start_year):(end_year),f1])/rowSums(pop.size.all[(start_year):(end_year),f1])) #overall prevalence in F, subpop 1
  prev.23.f <- mean((rowSums(inf.total[(start_year):(end_year),f2]) + rowSums(inf.total[(start_year):(end_year),f3])) / (rowSums(pop.size.all[(start_year):(end_year),f2])+ rowSums(pop.size.all[(start_year):(end_year),f3]))) # overall prevalence in F, subpops 2+3
  prev.y.1.f <- mean(rowSums(inf.total[(start_year):(end_year),y.f.1])/rowSums(pop.size.all[(start_year):(end_year),y.f.1])) #prevalence in F,subpop1, age cat 1
  prev.y.2.f<- mean(rowSums(inf.total[(start_year):(end_year),y.f.2])/rowSums(pop.size.all[(start_year):(end_year),y.f.2])) #prevalence in F,subpop2, age cat 1
  prev.y.23.f <- mean(rowSums(inf.total[(start_year):(end_year),y.f.23])/rowSums(pop.size.all[(start_year):(end_year),y.f.23])) #prevalence in F,subpops 2&3, age cat 1
  prev.o.1.f <- mean(rowSums(inf.total[(start_year):(end_year),o.f.1])/rowSums(pop.size.all[(start_year):(end_year),o.f.1])) #prevalence in F,subpop1, age cat 2
  prev.o.2.f<- mean(rowSums(inf.total[(start_year):(end_year),o.f.2])/rowSums(pop.size.all[(start_year):(end_year),o.f.2])) #prevalence in F,subpop2, age cat 2
  prev.o.23.f <- mean(rowSums(inf.total[(start_year):(end_year),o.f.23])/rowSums(pop.size.all[(start_year):(end_year),o.f.23])) #prevalence in F,subpops 2&3, age cat 1
  prev.y.1.m <- mean(rowSums(inf.total[(start_year):(end_year),y.m.1])/rowSums(pop.size.all[(start_year):(end_year),y.m.1])) #prevalence in M,subpop1, age cat 1
  prev.y.2.m<- mean(rowSums(inf.total[(start_year):(end_year),y.m.2])/rowSums(pop.size.all[(start_year):(end_year),y.m.2])) #prevalence in M,subpop2, age cat 1
  prev.y.3.m <- mean(rowSums(inf.total[(start_year):(end_year),y.m.3])/rowSums(pop.size.all[(start_year):(end_year),y.m.3])) #prevalence in M,subpops 2&3, age cat 1
  prev.o.1.m <- mean(rowSums(inf.total[(start_year):(end_year),o.m.1])/rowSums(pop.size.all[(start_year):(end_year),o.m.1])) #prevalence in M,subpop1, age cat 2
  prev.o.2.m<- mean(rowSums(inf.total[(start_year):(end_year),o.m.2])/rowSums(pop.size.all[(start_year):(end_year),o.m.2])) #prevalence in M,subpop2, age cat 2
  prev.o.3.m<- mean(rowSums(inf.total[(start_year):(end_year),o.m.3])/rowSums(pop.size.all[(start_year):(end_year),o.m.3])) #prevalence in M,subpop2, age cat 2
  prev.msw <- mean(rowSums(inf.total[(start_year):(end_year),msw])/rowSums(pop.size.all[(start_year):(end_year),msw])) #overall prevalence in msw
  # prev.msm <- sum(inf.total[(cal.start+cal.period+1),m4])/sum(pop.size.all[(cal.start+cal.period+1),m4]) #prevalence in msm in 2014
  prev.msm <- mean(rowSums(inf.total[(start_year):(end_year),m4]) / rowSums(pop.size.all[(start_year):(end_year), m4]))

  # Given the index of a population, named as y.m.1, y.f.1, o.m.1, etc.
  # calculate the prevalence of the described group.
  # If they're male, include the correct proportion of MSM.
  # These should evaluate to true:
  #   calculate_adjusted_prevalence(y.f.1) == mean(rowSums(inf.total[(start_year):(end_year),y.f.1])/(rowSums(pop.size.all[(start_year):(end_year),y.f.1])))
  #   calculate_adjusted_prevalence(y.m.1) == mean((rowSums(inf.total[(start_year):(end_year),y.m.1])+rowSums(inf.total[(start_year):(end_year), y.m.4])*p.s.1)/(rowSums(pop.size.all[(start_year):(end_year),y.m.1])+rowSums(pop.size.all[(start_year):(end_year), y.m.4])*p.s.1))
  calculate_adjusted_prevalence <- function(popindex) {

    numerator <- rowSums(inf.total[(start_year):(end_year), popindex])
    denominator <- rowSums(pop.size.all[(start_year):(end_year), popindex])

    popname <- deparse(substitute(popindex))
    popname_i <- gsub("[^0-9]", "", popname)

    if (grepl("m", popname)) {
      if(grepl("y", popname)) {
        age <- "y"
      } else if (grepl("o", popname)) {
        age <- "o"
      }
      msm_pop <- paste0(age, ".m.4")
      msm_pop <- eval(parse(text=msm_pop))
      msm_proportion <- paste0("p.s.", popname_i)
      msm_proportion <- eval(parse(text=msm_proportion))
      numerator <- numerator + rowSums(inf.total[(start_year):(end_year), msm_pop])*msm_proportion
      denominator <- denominator + rowSums(pop.size.all[(start_year):(end_year), msm_pop])*msm_proportion
    }
    return(mean(numerator/denominator))
  }

  # Prevalence rates with MSM Redistributed (for Male prevalence rates)
  prev.y.f.1.adj <- calculate_adjusted_prevalence(y.f.1)
  prev.y.m.1.adj <- calculate_adjusted_prevalence(y.m.1)
  prev.y.f.2.adj <- calculate_adjusted_prevalence(y.f.2)
  prev.y.m.2.adj <- calculate_adjusted_prevalence(y.m.2)
  prev.y.f.3.adj <- calculate_adjusted_prevalence(y.f.3)
  prev.y.m.3.adj <- calculate_adjusted_prevalence(y.m.3)

  prev.o.f.1.adj <- calculate_adjusted_prevalence(o.f.1)
  prev.o.m.1.adj <- calculate_adjusted_prevalence(o.m.1)
  prev.o.f.2.adj <- calculate_adjusted_prevalence(o.f.2)
  prev.o.m.2.adj <- calculate_adjusted_prevalence(o.m.2)
  prev.o.f.3.adj <- calculate_adjusted_prevalence(o.f.3)
  prev.o.m.3.adj <- calculate_adjusted_prevalence(o.m.3)

  # This ordering is meant to match exactly the ordering of the
  # data imported as gc_env$nhanes.updated.dat for use in the
  # prevalence likelihood function.
  msm_redistributed_prevalence_rates <- c(
    prev.y.f,
    prev.y.m,
    prev.y.f.1.adj,
    prev.y.m.1.adj,
    prev.y.f.2.adj,
    prev.y.m.2.adj,
    prev.y.f.3.adj,
    prev.y.m.3.adj,
    prev.o.f,
    prev.o.m,
    prev.o.f.1.adj,
    prev.o.m.1.adj,
    prev.o.f.2.adj,
    prev.o.m.2.adj,
    prev.o.f.3.adj,
    prev.o.m.3.adj
  )

  names(msm_redistributed_prevalence_rates) <- c(
    "prev.y.f",
    "prev.y.m",
    "prev.y.f.1.adj",
    "prev.y.m.1.adj",
    "prev.y.f.2.adj",
    "prev.y.m.2.adj",
    "prev.y.f.3.adj",
    "prev.y.m.3.adj",
    "prev.o.f",
    "prev.o.m",
    "prev.o.f.1.adj",
    "prev.o.m.1.adj",
    "prev.o.f.2.adj",
    "prev.o.m.2.adj",
    "prev.o.f.3.adj",
    "prev.o.m.3.adj"
  )


  ###calculate incidence for period covered by calibration
  inc <- annual_rates(sol[,1+inc.index])  #annual_rates assumes model is outputting ANNUAL values -- would need to change function if producing weekly/monthly outputs
  inc.s.a <- (sapply(1:((ncol(inc)-4)/2),function(x){rowSums(inc[,1:2+(x-1)*2])})) #incident CASES for given age, sex, and subpopulation group
  cal.start.diag <- cal.end - nrow(diag.age.sex.rate) + 1   #calculate years of data based on input (all data end in year 2016)
  n.inc <- as.vector(inc.s.a[(cal.start.diag):(cal.end),] ) #incident cases, 2000-2016
  ###calculate diagnosed cases by age, sex, and race/ethnicity, and redistribute cases in MSM to other male populations (based on proportionate sizes of these subpops)
  diagn <- annual_rates(sol[,1+diag.index]) #get annual number of reported cases (model outputs cumulative cases)
  diag.s.a <- (sapply(1:((ncol(diagn)-4)/2),function(x){rowSums(diagn[,1:2+(x-1)*2])})) #diagnosed CASES for given age, sex, and subpopulation group
  diag.s.a.adj <-diag.s.a[,1:6] + cbind(diag.s.a[,7:8]*p.s.1, diag.s.a[,7:8]*p.s.2, diag.s.a[,7:8]*p.s.3) #redistribute MSM cases among other male subpops
  diag.s.a.adj <- cbind(diag.s.a.adj,diag.s.a[,9:14]) #append female cases to adjusted male cases
  diag.s.a.adj.rate <-t(t(diag.s.a.adj)/pop.sub.adj) #diagnosed case rate, with MSM cases redistributed among other male subpops
  diag.s.a.msm.rate <- t(t(diag.s.a[,7:8])/pop.msm) #diagnosed case rate in MSM


  # we're going to calibrate to the diag.subpop.rate
  # the diag.subpop.rate is going to get unlisted as columns in order
  # This is the flattened data we're going to calibrate to:
  # diag.subpop.rate_flat <- unlist(c(diag.subpop.rate))
  # names(diag.subpop.rate_flat) <- rep(names(diag.subpop.rate), each = 7)

  #reported case rates by age and sex
  n.diag.age.sex.time <- c(rowSums(diagn[((cal.start.diag):(cal.end)),y.m])/sum(n.i[y.m]),
                           rowSums(diagn[((cal.start.diag):(cal.end)),o.m])/sum(n.i[o.m]),
                           rowSums(diagn[((cal.start.diag):(cal.end)),y.f])/sum(n.i[y.f]),
                           rowSums(diagn[((cal.start.diag):(cal.end)),o.f])/sum(n.i[o.f]))

  #### calculate relative reported case rate for time period for which have data by race/ethnicity, using overall rates by age and sex as comparator
  cal.start.rr <-  cal.end - nrow(diag.rr)/8 + 1   #calculate years of data on subpop RR based on input (all data end in year 2016)
  length.rr <- nrow(diag.rr)/8 #number of data points for subpop RR data

  #reported case rates by age and sex, for time period for which have race/ethnicity RR data
  #using this as denominator for rr diag in black and Hispanic populations
  n.diag.rr <- c(rowSums(diagn[((cal.start.rr):(cal.end)),y.m])/sum(n.i[y.m]),
                           rowSums(diagn[((cal.start.rr):(cal.end)),o.m])/sum(n.i[o.m]),
                           rowSums(diagn[((cal.start.rr):(cal.end)),y.f])/sum(n.i[y.f]),
                           rowSums(diagn[((cal.start.rr):(cal.end)),o.f])/sum(n.i[o.f]))

  rr.diag.y.m.1 <- diag.s.a.adj.rate[(cal.start.rr):(cal.end),1]/ n.diag.rr[1:length.rr]
  rr.diag.o.m.1 <- diag.s.a.adj.rate[(cal.start.rr):(cal.end),2]/ n.diag.rr[(length.rr+1):(2*length.rr)]
  rr.diag.y.m.3 <- diag.s.a.adj.rate[(cal.start.rr):(cal.end),5]/n.diag.rr[1:length.rr]
  rr.diag.o.m.3 <- diag.s.a.adj.rate[(cal.start.rr):(cal.end),6]/ n.diag.rr[(length.rr+1):(2*length.rr)]
  rr.diag.y.f.1 <- diag.s.a.adj.rate[(cal.start.rr):(cal.end),7]/ n.diag.rr[(2*length.rr+1):(3*length.rr)]
  rr.diag.o.f.1 <- diag.s.a.adj.rate[(cal.start.rr):(cal.end),8]/  n.diag.rr[(3*length.rr+1):(4*length.rr)]
  rr.diag.y.f.3 <- diag.s.a.adj.rate[(cal.start.rr):(cal.end),11]/ n.diag.rr[(2*length.rr+1):(3*length.rr)]
  rr.diag.o.f.3 <- diag.s.a.adj.rate[(cal.start.rr):(cal.end),12]/ n.diag.rr[(3*length.rr+1):(4*length.rr)]

  #calculate proportion of diagnosed cases with symptoms
  years.symp <- 7 #years of data averaging for pSymptomatic data, currently using 2010-2016
  cal.start.symp <- cal.end +1 - years.symp
  symp.diag <- annual_rates(sol[,(1+symp.index)])/diagn #proportion diagnosed cases that are symptomatic
  cal.end.symp <- cal.end
  if (gc_env$site == 'BA') cal.end.symp <- cal.end.symp + 1
  p.symp.diag <- c(weighted.mean(symp.diag[(cal.start.symp):(cal.end.symp), y.m.msw]),
                   weighted.mean(symp.diag[(cal.start.symp):(cal.end.symp), o.m.msw]),
                   weighted.mean(symp.diag[(cal.start.symp):(cal.end.symp), y.m.4]),
                   weighted.mean(symp.diag[(cal.start.symp):(cal.end.symp) ,y.m.msw]),
                   weighted.mean(symp.diag[(cal.start.symp):(cal.end.symp) ,y.f]),
                   weighted.mean(symp.diag[(cal.start.symp):(cal.end.symp) ,o.f]))

  #calculate proportion of male cases who are MSM
  # since we have new SSuN data from baltimore, we should include 2017 in the
  # p.diag.msm calculation for Baltimore but not for SF
  diag.msm.cal.end <- cal.end
  if (gc_env$site == 'BA') {
    diag.msm.cal.end <- cal.start + 15 # 2017
  } else diag.msm.cal.end <- cal.start + 14 # 2016
  # cal.start.msm <- cal.end +1 - nrow(p.msm.dat)/2
  cal.start.msm <- cal.start + 8 # 2010
  p.diag.msm <-  c(rowSums(diagn[(cal.start.msm):(diag.msm.cal.end),y.m.4]),rowSums(diagn[(cal.start.msm):(diag.msm.cal.end),o.m.4]) ) /
                 c(rowSums(diagn[(cal.start.msm):(diag.msm.cal.end),y.m]),rowSums(diagn[(cal.start.msm):(diag.msm.cal.end),o.m]) )

  rate.diag.subpop.time <- as.vector(diag.s.a.adj.rate[(cal.start.rr):(cal.end),] )# rate of diagnosed cases by subpop, sex, age
  rate.diag.subpop.time.msm<-as.vector(diag.s.a.msm.rate[(cal.start.rr):(cal.end),]) #rate of diagnosed cases, MSM

  fit.prev.extra <-
    c(
      prev.m.adj = prev.m.adj,
      prev.f = prev.f,
      prev.y.m.adj = prev.y.m.adj,
      prev.y.f = prev.y.f,
      prev.1.adj = prev.1.adj,
      prev.23.adj = prev.23.adj,
      prev.age.subpop.m.adj.msm = prev.age.subpop.m.adj.msm[c(1, 3, 2, 4)],
      prev.y.1.f = prev.y.1.f,
      prev.y.23.f = prev.y.23.f,
      prev.o.1.f = prev.o.1.f,
      prev.o.23.f = prev.o.23.f,
      prev.msm = prev.msm
    )
  fit.prev.extra.no.msm <-
    c(
      prev.msw = prev.msw,
      prev.f = prev.f,
      prev.y.msw = prev.y.msw,
      prev.y.f = prev.y.f,
      prev.1.no.msm = prev.1.no.msm,
      prev.23.no.msm = prev.23.no.msm,
      prev.y.1.m = prev.y.1.m,
      prev.y.2.m = prev.y.2.m,
      prev.y.3.m = prev.y.3.m,
      prev.o.1.m = prev.o.1.m,
      prev.o.2.m = prev.o.2.m,
      prev.o.3.m = prev.o.3.m,
      prev.y.1.f = prev.y.1.f,
      prev.y.2.f = prev.y.2.f,
      prev.y.23.f = prev.y.23.f,
      prev.o.1.f = prev.o.1.f,
      prev.o.2.f = prev.o.2.f,
      prev.o.f = prev.o.f,
      prev.msm = prev.msm
    )
  fit.diag.rr <-
    c(
       rr.diag.y.m.1,
       rr.diag.o.m.1,
       rr.diag.y.m.3,
       rr.diag.o.m.3,
       rr.diag.y.f.1,
       rr.diag.o.f.1,
       rr.diag.y.f.3,
       rr.diag.o.f.3
    )
  fit.symp <-p.symp.diag
  list(
    fit.symp = fit.symp,
    inf.subpop.age = inf.subpop.age[1:14],
    fit.prev.extra = fit.prev.extra,
    fit.prev.extra.no.msm = fit.prev.extra.no.msm,
    fit.diag.rate = n.diag.age.sex.time,
    fit.diag.subpop = rate.diag.subpop.time,
    fit.diag.msm = rate.diag.subpop.time.msm ,
    p.diag.msm = p.diag.msm,
    inc = n.inc,
    diag.rr = fit.diag.rr,
    msm_redistributed_prevalence_rates = msm_redistributed_prevalence_rates
  )
}
