#########################################################
### functions for calculating time-varying parameters ###
#########################################################

#' Import the Bezier function from the bezier package, and not Hmisc
bezier <- bezier::bezier

# bezier function for screening rate/reporting, constrained to be between 0 and
# 1 function takes 4 parameters: 1 = start point (value at t=0), 2=end point
# (value at t=end), 3&4=internal control points produces annual estimates of
# screening/reporting rate for length specified (the calibration period)
bezier_fun<- function(screen.bezier){
  cal.period <- gc_env$cal.period
  bez.a<-screen.bezier[1]
  bez.d<-screen.bezier[2]
  bez.b <- screen.bezier[3]
  bez.c <- screen.bezier[4]
  y0=bez.a
  y1=bez.b
  y2=bez.c
  y3=bez.d
  p0 =  y0
  p1 = ( -5*y0 + 18*y1 -  9*y2 + 2*y3) / 6
  p2 = (  2*y0 -  9*y1 + 18*y2 - 5*y3) / 6
  p3 = y3
  if (exists('gc_testing', envir = .GlobalEnv) && 'first_three_bezier_equal' %in% gc_testing) {
    p2 <- p1 <- p0
  }
  p<-c(p0,p1,p2,p3)
  t<-seq(0,1, length.out=(cal.period))
  #this is the function from the bezier package - Hmisc also has a bezier
  #function that is NOT the one to use
  bezier::bezier(t,p)
}

#calculate control points b and c for bezier curves describing screening and reporting
#see model technical appendix for additonal details
update_ctrl <- function(theta) {

  #internal control points for screening rate in Hispanic and other females, 15-24y
  screen.f1.b = ilogit(theta["logit.screen.f1.a"]) + (ilogit(theta["logit.screen.f1.d"]) -
    ilogit(theta["logit.screen.f1.a"])) * ilogit(theta["logit.rand.screen.f1.b"])
  screen.f1.c = screen.f1.b + (ilogit(theta["logit.screen.f1.d"]) - screen.f1.b) *
    ilogit(theta["logit.rand.screen.f1.c"])
  #internal control points for screening rate in Hispanic and other females, 25-39y
  screen.f2.b = ilogit(theta["logit.screen.f2.a"]) + (ilogit(theta["logit.screen.f2.d"]) -
    ilogit(theta["logit.screen.f2.a"])) * ilogit(theta["logit.rand.screen.f2.b"])
  screen.f2.c = screen.f2.b + (ilogit(theta["logit.screen.f2.d"]) - screen.f2.b) *
    ilogit(theta["logit.rand.screen.f2.c"])

  #internal control points for screening rate in black females, 15-24y
  screen.f1.1.b = ilogit(theta["logit.screen.f1.1.a"]) + (ilogit(theta["logit.screen.f1.1.d"]) -
    ilogit(theta["logit.screen.f1.1.a"])) * ilogit(theta["logit.rand.screen.f1.1.b"])
  screen.f1.1.c = screen.f1.1.b + (ilogit(theta["logit.screen.f1.1.d"]) -
    screen.f1.1.b) * ilogit(theta["logit.rand.screen.f1.1.c"])
  #internal control points for screening rate in black females, 25-39y
  screen.f2.1.b = ilogit(theta["logit.screen.f2.1.a"]) + (ilogit(theta["logit.screen.f2.1.d"]) -
    ilogit(theta["logit.screen.f2.1.a"])) * ilogit(theta["logit.rand.screen.f2.1.b"])
  screen.f2.1.c = screen.f2.1.b + (ilogit(theta["logit.screen.f2.1.d"]) -
                                     screen.f2.1.b) * ilogit(theta["logit.rand.screen.f2.1.c"])

  #internal control points for screening rate in Hispanic and other males, 15-24y
  screen.m1.b = ilogit(theta["logit.screen.m1.a"]) + (ilogit(theta["logit.screen.m1.d"]) -
        ilogit(theta["logit.screen.m1.a"])) * ilogit(theta["logit.rand.screen.m1.b"])
  screen.m1.c = screen.m1.b + (ilogit(theta["logit.screen.m1.d"]) - screen.m1.b) *
    ilogit(theta["logit.rand.screen.m1.c"])
  #internal control points for screening rate in Hispanic and other males, 25-39y
  screen.m2.b = ilogit(theta["logit.screen.m2.a"]) + (ilogit(theta["logit.screen.m2.d"]) -
                ilogit(theta["logit.screen.m2.a"])) * ilogit(theta["logit.rand.screen.m2.b"])
  screen.m2.c = screen.m2.b + (ilogit(theta["logit.screen.m2.d"]) - screen.m2.b) *
    ilogit(theta["logit.rand.screen.m2.c"])

  #internal control points for screening rate in black males, 15-24y
  screen.m1.1.b = ilogit(theta["logit.screen.m1.1.a"]) + (ilogit(theta["logit.screen.m1.1.d"]) -
                  ilogit(theta["logit.screen.m1.1.a"])) * ilogit(theta["logit.rand.screen.m1.1.b"])
  screen.m1.1.c = screen.m1.1.b + (ilogit(theta["logit.screen.m1.1.d"]) -
                  screen.m1.1.b) * ilogit(theta["logit.rand.screen.m1.1.c"])
  #internal control points for screening rate in black males, 25-39y
  screen.m2.1.b = ilogit(theta["logit.screen.m2.1.a"]) + (ilogit(theta["logit.screen.m2.1.d"]) -
                  ilogit(theta["logit.screen.m2.1.a"])) * ilogit(theta["logit.rand.screen.m2.1.b"])
  screen.m2.1.c = screen.m2.1.b + (ilogit(theta["logit.screen.m2.1.d"]) -
                  screen.m2.1.b) * ilogit(theta["logit.rand.screen.m2.1.c"])

  #internal control points for screening rate in MSM, 15-24y
  screen.msm1.b = ilogit(theta["logit.screen.msm1.a"]) + (ilogit(theta["logit.screen.msm1.d"]) -
                  ilogit(theta["logit.screen.msm1.a"])) * ilogit(theta["logit.rand.screen.msm1.b"])
  screen.msm1.c = screen.msm1.b + (ilogit(theta["logit.screen.msm1.d"]) -
                  screen.msm1.b) * ilogit(theta["logit.rand.screen.msm1.c"])
  #internal control points for screening rate in MSM, 25-39y
  screen.msm2.b = ilogit(theta["logit.screen.msm2.a"]) + (ilogit(theta["logit.screen.msm2.d"]) -
                  ilogit(theta["logit.screen.msm2.a"])) * ilogit(theta["logit.rand.screen.msm2.b"])
  screen.msm2.c = screen.msm2.b + (ilogit(theta["logit.screen.msm2.d"]) -
                  screen.msm2.b) * ilogit(theta["logit.rand.screen.msm2.c"])

  #save internal control points for screening
  screen.bez.ctrl <-
    matrix(
      c(
        screen.m1.1.b,
        screen.m1.1.c,
        screen.m1.b,
        screen.m1.c,
        screen.m2.1.b,
        screen.m2.1.c,
        screen.m2.b,
        screen.m2.c,
        screen.msm1.b,
        screen.msm1.c,
        screen.msm2.b,
        screen.msm2.c,
        screen.f1.1.b,
        screen.f1.1.c,
        screen.f1.b,
        screen.f1.c,
        screen.f2.1.b,
        screen.f2.1.c,
        screen.f2.b,
        screen.f2.c
      ),
      ncol = 2,
      byrow = T
    )
  rownames(screen.bez.ctrl)<-c("m1.1","m1","m2.1","m2","msm1","msm2","f1.1","f1","f2.1","f2")
  colnames(screen.bez.ctrl)<-c("b","c")

  #internal control points for reporting parameter
  rep.symp.b = ilogit(theta["logit.rep.symp.a"]) + (ilogit(theta["logit.rep.symp.d"]) -
               ilogit(theta["logit.rep.symp.a"])) * ilogit(theta["logit.rand.rep.symp.b"])
  rep.symp.c = rep.symp.b + (ilogit(theta["logit.rep.symp.d"]) -
               rep.symp.b) * ilogit(theta["logit.rand.rep.symp.c"])
  rep.bez.bc <- c(rep.symp.b, rep.symp.c)
  gc_env$rep.bez.bc <- rep.bez.bc

  return(list(screen.bez.ctrl, rep.bez.bc))

}

#calculate control points b and c for bezier curves for plotting model inputs
prior_ctrl <- function(bez) {
  bez.b <- bez[1] + (bez[2]-bez[1])*bez[3]
  bez.c <- bez.b + (bez[2]-bez.b)*bez[4]
  prior <- unlist(c(bez[1], bez[2], b=bez.b, c=bez.c) )
  return(prior)
}

#calculate screening rate from base screening rate x rr screening in high sexual
#activity group, for each row of matrix
screen_fun <-function(screen.param,rr.screen) {
  x=NULL
  for (i in 1:length(screen.param)) {
    x<-c(x,screen.param[i], screen.param[i]*rr.screen)
  }
  return(c(x,rep(0,4))) #add on 0's for the empty i=4 group for females
}


#calculate time-varying tranmission relative risk for MSM
behav_fun<- function(behav.m){
  cal.period <- gc_env$cal.period
  t<-seq(from=0,to=cal.period, by=1)
  t <- t / cal.period
  return(behav.m*t)
}
