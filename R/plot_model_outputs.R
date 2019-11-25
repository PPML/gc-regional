###########################################################################
### visualize prior and posterior distributions and calibration targets ###
###########################################################################

#run this after calibration, with burned, trimmed, and merged trace


#' @import ggplot2
get_legend<-function(myggplot){   # function to get legend from plot
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#' @title Plot Calibration Targets and Priors/Posteriors
#'
#' @export
#' @import reshape2
#' @import grid
#' @import gridExtra
#' @import ggplot2
plot_posteriors <- function(output_directory, filename = paste0('Calibration_Plots_', gc_env$site, ".pdf")) {
  require(reshape2)
  # Retrieve all environment variables from gc_env
  for (varname in ls(envir=gc_env)) {
    assign(varname, gc_env[[varname]])
  }

  sitename <- if(site=="SF") "San Francisco" else "Baltimore"

  blankPlot <- ggplot()+geom_blank(aes(1,1)) +  # make a blank plot to use as a placeholder where necessary
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank()
    )

  #prep input data
  diag.age.sex.rate <- if(!"year" %in% colnames(as.data.frame(diag.age.sex.rate))) tibble::rownames_to_column(as.data.frame(diag.age.sex.rate), "year")  else diag.age.sex.rate
  diag.dat <- melt(diag.age.sex.rate[,1:5],  id.vars = "year", measured.vars=c("m_y", "m_o", "f_y","f_o"))
  diag.dat$cat <- ifelse(diag.dat$variable=="m_y", "M 15-24 y", ifelse(diag.dat$variable=="m_o", "M 25-39 y", ifelse(diag.dat$variable=="f_y", "F 15-24 y", "F 25-39 y")))
  p.msm.dat$cat <- ifelse(p.msm.dat$age_cat==1, "15-24 y", "25-39 y" )
  diag.subpop.rate <- if(!"year" %in% colnames(diag.subpop.rate)) tibble::rownames_to_column(diag.subpop.rate, "year")  else diag.subpop.rate
  subpop.dat <- melt(diag.subpop.rate,  id.vars = "year", measured.vars=list[2:ncol(diag.subpop.rate)])
  subpop.dat$cat <- ifelse(subpop.dat$variable=="y.m.1", "Black M 15-24 y",
  ifelse(subpop.dat$variable=="y.m.2", "Other M 15-24 y",
  ifelse(subpop.dat$variable=="y.m.3", "Hispanic M 15-24 y",
  ifelse(subpop.dat$variable=="o.m.1", "Black M 25-39 y",
  ifelse(subpop.dat$variable=="o.m.2", "Other M 25-39 y",
  ifelse(subpop.dat$variable=="o.m.3", "Hispanic M 25-39 y",
  ifelse(subpop.dat$variable=="y.f.1", "Black F 15-24 y",
  ifelse(subpop.dat$variable=="y.f.2", "Other F 15-24 y",
  ifelse(subpop.dat$variable=="y.f.3", "Hispanic F 15-24 y",
  ifelse(subpop.dat$variable=="o.f.1", "Black F 25-39 y",
  ifelse(subpop.dat$variable=="o.f.2", "Other F 25-39 y", "Hispanic F 25-39 y")))))))))))

	subpop.dat$sex <- ifelse(subpop.dat$variable=="y.m.1"|subpop.dat$variable=="y.m.2"|subpop.dat$variable=="y.m.3"|
															 subpop.dat$variable=="o.m.1"|subpop.dat$variable=="o.m.2"|subpop.dat$variable=="o.m.3", "M", "F")


  subpop.dat$age <- ifelse(subpop.dat$variable=="y.m.1"|subpop.dat$variable=="y.m.2"|subpop.dat$variable=="y.m.3"|
                             subpop.dat$variable=="y.f.1"|subpop.dat$variable=="y.f.2"|subpop.dat$variable=="y.f.3", "15-24 y", "25-39 y")

  subpop.dat$Population <- ifelse(subpop.dat$variable=="y.m.1"|subpop.dat$variable=="o.m.1"|subpop.dat$variable=="y.f.1"|subpop.dat$variable=="o.f.1", "Black",
                                  ifelse(subpop.dat$variable=="y.m.2"|subpop.dat$variable=="o.m.2"|subpop.dat$variable=="y.f.2"|subpop.dat$variable=="o.f.2", "Other", "Hispanic"))

  # prep outputs for plotting
  # here y=15-24y, o=25-39y, m=male, f=female, 1=black, 2=other, 3=Hispanic, msm=men who have sex with men
  length.rr <- nrow(diag.rr)/8
  # out.prev <- melt(as.matrix(subset(pred, select=prev.fit.prev.extra1:prev.fit.prev.extra15)))  #model prevalence with extra categories, with MSM in calculations
  out.prev <- pred[, grep('prev.msm_redistributed_prevalence_rates|prev.fit.prev.extra.prev.msm', names(pred), value=T)]
  out.diag.y.m.1 <- melt(as.matrix(subset(pred, select=prev.fit.diag.subpop1:(prev.fit.diag.subpop1+length.rr-1)))) #reported cases by subpop
  out.diag.o.m.1 <- melt(as.matrix(subset(pred, select=(prev.fit.diag.subpop1+length.rr):(prev.fit.diag.subpop1+2*length.rr-1))))
  out.diag.y.m.2 <- melt(as.matrix(subset(pred, select=(prev.fit.diag.subpop1+2*length.rr):(prev.fit.diag.subpop1+3*length.rr-1))))
  out.diag.o.m.2 <- melt(as.matrix(subset(pred, select=(prev.fit.diag.subpop1+3*length.rr):(prev.fit.diag.subpop1+4*length.rr-1))))
  out.diag.y.m.3 <- melt(as.matrix(subset(pred, select=(prev.fit.diag.subpop1+4*length.rr):(prev.fit.diag.subpop1+5*length.rr-1))))
  out.diag.o.m.3 <- melt(as.matrix(subset(pred, select=(prev.fit.diag.subpop1+5*length.rr):(prev.fit.diag.subpop1+6*length.rr-1))))
  out.diag.y.f.1 <- melt(as.matrix(subset(pred, select=(prev.fit.diag.subpop1+6*length.rr):(prev.fit.diag.subpop1+7*length.rr-1))))
  out.diag.o.f.1 <- melt(as.matrix(subset(pred, select=(prev.fit.diag.subpop1+7*length.rr):(prev.fit.diag.subpop1+8*length.rr-1))))
  out.diag.y.f.2 <- melt(as.matrix(subset(pred, select=(prev.fit.diag.subpop1+8*length.rr):(prev.fit.diag.subpop1+9*length.rr-1))))
  out.diag.o.f.2 <- melt(as.matrix(subset(pred, select=(prev.fit.diag.subpop1+9*length.rr):(prev.fit.diag.subpop1+10*length.rr-1))))
  out.diag.y.f.3 <- melt(as.matrix(subset(pred, select=(prev.fit.diag.subpop1+10*length.rr):(prev.fit.diag.subpop1+11*length.rr-1))))
  out.diag.o.f.3 <- melt(as.matrix(subset(pred, select=(prev.fit.diag.subpop1+11*length.rr):(prev.fit.diag.subpop1+12*length.rr-1))))
  out.diag.y.msm <- melt(as.matrix(subset(pred, select=prev.fit.diag.msm1:(prev.fit.diag.msm1+length.rr-1))))
  out.diag.o.msm <- melt(as.matrix(subset(pred, select=(prev.fit.diag.msm1+length.rr):(prev.fit.diag.msm1+2*length.rr-1))))

  length.diag <- nrow(diag.age.sex.rate)
  out.diag.y.m <-melt(as.matrix(subset(pred, select=prev.fit.diag.rate1:(prev.fit.diag.rate1+length.diag-1)))) # reported cases rates by age and sex only
  out.diag.o.m <-melt(as.matrix(subset(pred, select=(prev.fit.diag.rate1+length.diag):(prev.fit.diag.rate1+2*length.diag-1))))
  out.diag.y.f <-melt(as.matrix(subset(pred, select=(prev.fit.diag.rate1+2*length.diag):(prev.fit.diag.rate1+3*length.diag-1))))
  out.diag.o.f<-melt(as.matrix(subset(pred, select=(prev.fit.diag.rate1+3*length.diag):(prev.fit.diag.rate1+4*length.diag-1))))

  out.inc.y.m.1 <-melt(100*as.matrix(subset(pred, select=prev.inc1:(prev.inc1+length.diag-1)))/sum(n.i[y.m.1])) #model incidence by age, sex, and subpopulation
  out.inc.o.m.1 <-melt(100*as.matrix(subset(pred, select=(prev.inc1+length.diag):(prev.inc1+2*length.diag-1)))/sum(n.i[o.m.1]))
  out.inc.y.m.2 <-melt(100*as.matrix(subset(pred, select=(prev.inc1+2*length.diag):(prev.inc1+3*length.diag-1)))/sum(n.i[y.m.2]))
  out.inc.o.m.2 <-melt(100*as.matrix(subset(pred, select=(prev.inc1+3*length.diag):(prev.inc1+4*length.diag-1)))/sum(n.i[o.m.2]))
  out.inc.y.m.3 <-melt(100*as.matrix(subset(pred, select=(prev.inc1+4*length.diag):(prev.inc1+5*length.diag-1)))/sum(n.i[y.m.3]))
  out.inc.o.m.3 <-melt(100*as.matrix(subset(pred, select=(prev.inc1+5*length.diag):(prev.inc1+6*length.diag-1)))/sum(n.i[o.m.3]))
  out.inc.y.msm <-melt(100*as.matrix(subset(pred, select=(prev.inc1+6*length.diag):(prev.inc1+7*length.diag-1)))/sum(n.i[y.m.4]))
  out.inc.o.msm <-melt(100*as.matrix(subset(pred, select=(prev.inc1+7*length.diag):(prev.inc1+8*length.diag-1)))/sum(n.i[o.m.4]))
  out.inc.y.f.1 <-melt(100*as.matrix(subset(pred, select=(prev.inc1+8*length.diag):(prev.inc1+9*length.diag-1)))/sum(n.i[y.f.1]))
  out.inc.o.f.1 <-melt(100*as.matrix(subset(pred, select=(prev.inc1+9*length.diag):(prev.inc1+10*length.diag-1)))/sum(n.i[o.f.1]))
  out.inc.y.f.2 <-melt(100*as.matrix(subset(pred, select=(prev.inc1+10*length.diag):(prev.inc1+11*length.diag-1)))/sum(n.i[y.f.2]))
  out.inc.o.f.2 <-melt(100*as.matrix(subset(pred, select=(prev.inc1+11*length.diag):(prev.inc1+12*length.diag-1)))/sum(n.i[o.f.2]))
  out.inc.y.f.3 <-melt(100*as.matrix(subset(pred, select=(prev.inc1+12*length.diag):(prev.inc1+13*length.diag-1)))/sum(n.i[y.f.3]))
  out.inc.o.f.3 <-melt(100*as.matrix(subset(pred, select=(prev.inc1+13*length.diag):(prev.inc1+14*length.diag-1)))/sum(n.i[o.f.3]))

  out.rr.diag.m.y.1 <- melt(as.matrix(subset(pred, select=prev.diag.rr1:(prev.diag.rr1+length.rr-1)))) #relative risk of being a reported case by subpop, age, sex (here 1=black, 2=Hispanic)
  out.rr.diag.m.o.1 <- melt(as.matrix(subset(pred, select=(prev.diag.rr1+length.rr):(prev.diag.rr1+2*length.rr-1))))
  out.rr.diag.m.y.2 <- melt(as.matrix(subset(pred, select=(prev.diag.rr1+2*length.rr):(prev.diag.rr1+3*length.rr-1))))
  out.rr.diag.m.o.2 <- melt(as.matrix(subset(pred, select=(prev.diag.rr1+3*length.rr):(prev.diag.rr1+4*length.rr-1))))
  out.rr.diag.f.y.1 <- melt(as.matrix(subset(pred, select=(prev.diag.rr1+4*length.rr):(prev.diag.rr1+5*length.rr-1))))
  out.rr.diag.f.o.1 <- melt(as.matrix(subset(pred, select=(prev.diag.rr1+5*length.rr):(prev.diag.rr1+6*length.rr-1))))
  out.rr.diag.f.y.2 <- melt(as.matrix(subset(pred, select=(prev.diag.rr1+6*length.rr):(prev.diag.rr1+7*length.rr-1))))
  out.rr.diag.f.o.2 <- melt(as.matrix(subset(pred, select=(prev.diag.rr1+7*length.rr):(prev.diag.rr1+8*length.rr-1))))

  out.subpop <- melt(as.matrix(subset(pred, select=pred.s.dist1:pred.s.dist6))) #model subpopulation assortativity
  out.age <- melt(as.matrix(subset(pred, select=age.dist.all1:age.dist.all4))) #model age assortativity
  out.symp <- melt(as.matrix(subset(pred, select=prev.fit.symp1:prev.fit.symp6))) #model proportion of diagnosed cases that are symptomatic
  length.msm <- nrow(p.msm.dat)/2
  out.p.msm.y <- melt(as.matrix(subset(pred, select=prev.p.diag.msm1:(prev.p.diag.msm1+length.msm-1))))
  out.p.msm.o <- melt(as.matrix(subset(pred, select=(prev.p.diag.msm1+length.msm):(prev.p.diag.msm1+2*length.msm-1))))

  diag.data <- cbind(year=(seq(as.numeric(diag.age.sex.rate[1,"year"]),as.numeric(diag.age.sex.rate[nrow(diag.age.sex.rate), "year"]),1)), as.data.frame(diag.rate), cat=rep(c("y.m", "o.m", "y.f", "o.f"),each=nrow(diag.age.sex.rate) ))
  age.data <- as.data.frame(age.dist.dat[1:4,])
  p.msm.plot <- p.msm.ssun

  ### plot input data and model outputs ###

  # plot.prev <- ggplot(data=out.prev)+
  #   geom_pointrange(data=prev.extra.cat.dat,aes(x=(1:nrow(prev.extra.cat.dat)+0), y=prev, ymin=LCL, ymax=UCL), color="dodgerblue", shape=15, size=1, alpha=.3) +
  #   geom_boxplot(aes(x=Var2, y=value, fill=Var2), outlier.color="darkgrey",position="dodge") +
  #   theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1,size=8),axis.text.y=element_text(size=8), title=element_text(size=9)) +
  #   labs(title="Gonorrhea prevalence\n(not used in calibration)", x="\nPopulation", y="Prevalence (%)\n") +
  #   scale_x_discrete(labels=c("M all","F all", "M y", "F y", "Black MF all", "Hisp+other MF all","Black M y","Hisp+other M y","Black M o","Hisp+other M o", "Black F y", "Hisp+other F y", "Black F o", "Hisp+Other F o", "All MSM"))
  nhanes.updated.dat <- gc_env$nhanes.updated.dat

  shortname <- c("Y F", "Y M", "Y F Black", "Y M Black", "Y F Other", "Y M Other", "Y F Hisp", "Y M Hisp", "O F", "O M", "O F Black", "O M Black", "O F Other", "O M Other", "O F Hisp", "O M Hisp")
  # colnames(out.prev) <- gsub("prev.msm_redistributed_prevalence_rates.prev.", "", colnames(out.prev))
  nhanes_cols <- grepl("prev.msm_redistributed_prevalence_rates.prev", colnames(out.prev))
  colnames(out.prev)[nhanes_cols] <- shortname
  colnames(out.prev) <- gsub("prev.fit.prev.extra.prev.msm", "All MSM", colnames(out.prev))
  nhanes.updated.dat[['shortname']] <- shortname

  out.prev <- reshape2::melt(out.prev)

  plot.prev <- ggplot(data = out.prev) +
    geom_pointrange(data = nhanes.updated.dat,
                    aes(
                      x = shortname,
                      y = prev,
                      ymin = prev - prev_std,
                      ymax = prev + prev_std
                    ), color = 'red', shape = 15, alpha = 0.7) +
    geom_boxplot(aes(x=variable, y=value, fill = variable), width = 0.75, outlier.alpha=0.3, position = "dodge", alpha=.6) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=90, hjust=1, size=8), legend.position = 'none', axis.text.y=element_text(size=8), title=element_text(size=8)) +
    labs(title="Gonorrhea Prevalence", x="\nPopulation", y="Prevalence (%)\n")


  if (gc_env$site == "SF") {
    plot.prev <- plot.prev +
      geom_pointrange(data = data.frame(x="All MSM", y=.06), aes(x=x,y=y,ymax=y+0.012, ymin=y-0.012),  color="dodgerblue", shape=15, alpha=.7)
  }

  #get max values of reported cases for setting y-axis
  max.y.m <- 100* max(out.diag.y.m.1$value,out.diag.y.m.2$value,out.diag.y.m.3$value,diag.subpop.rate[,"y.m.1"],diag.subpop.rate[,"y.m.2"], diag.subpop.rate[,"y.m.3"], na.rm=TRUE)
  max.o.m <- 100*max(out.diag.o.m.1$value,out.diag.o.m.2$value,out.diag.o.m.3$value,diag.subpop.rate[,"o.m.1"], diag.subpop.rate[,"o.m.2"],diag.subpop.rate[,"o.m.3"], na.rm=TRUE)
  max.y.f <- 100*max(out.diag.y.f.1$value,out.diag.y.f.2$value,out.diag.y.f.3$value,diag.subpop.rate[,"y.f.1"], diag.subpop.rate[,"y.f.2"],diag.subpop.rate[,"y.f.3"], na.rm=TRUE)
  max.o.f <- 100*max(out.diag.o.f.1$value,out.diag.o.f.2$value,out.diag.o.f.3$value, diag.subpop.rate[,"o.f.1"],diag.subpop.rate[,"o.f.2"], diag.subpop.rate[,"o.f.3"], na.rm=TRUE)
  max.diag <- max(max.y.m, max.o.m, max.y.f, max.o.f) +0.1
  max.msm <- max(out.diag.y.msm$value,out.diag.o.msm$value, na.rm=TRUE)*100+0.1

  ### plot diagnosed cases by age and sex ###
  plot.diag.y.m <- ggplot(data=out.diag.y.m)+
    geom_line(aes(x=Var2, y=value*100, group=Var1), color="grey") +
    stat_summary(aes(x=Var2, y=value*100, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
    geom_point(data=diag.data[1:length.diag,],aes(x=1:length.diag, y=diag.rate*100), color="red", shape=15, size=2) +
    labs(title="All M: 15-24y", x="Year", y="Reported diagnoses per 100") +
    coord_cartesian(ylim=c(0,max.y.m))+
    scale_x_discrete(labels=seq((end.year-length.diag+1),end.year,1))

  plot.diag.o.m <- ggplot(data=out.diag.o.m)+
    geom_line(aes(x=Var2, y=value*100, group=Var1), color="grey") +
    stat_summary(aes(x=Var2, y=value*100, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
    geom_point(data=diag.data[(length.diag+1):(2*length.diag),],aes(x=1:length.diag, y=diag.rate*100), color="red", shape=15, size=2) +
    labs(title="All M: 25-39y", x="Year", y="Reported diagnoses per 100") +
    coord_cartesian(ylim=c(0,max.o.m))+
    scale_x_discrete(labels=seq((end.year-length.diag+1),end.year,1))

  plot.diag.y.f <- ggplot(data=out.diag.y.f)+
    geom_line(aes(x=Var2, y=value*100, group=Var1), color="grey") +
    stat_summary(aes(x=Var2, y=value*100, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
    geom_point(data=diag.data[(2*length.diag+1):(3*length.diag),],aes(x=1:length.diag, y=diag.rate*100), color="red", shape=15, size=2) +
    labs(title="All F: 15-24y", x="Year", y="Reported diagnoses per 100") +
    coord_cartesian(ylim=c(0,max.y.f))+
    scale_x_discrete(labels=seq((end.year-length.diag+1),end.year,1))

  plot.diag.o.f <- ggplot(data=out.diag.o.f)+
    geom_line(aes(x=Var2, y=value*100, group=Var1), color="grey") +
    stat_summary(aes(x=Var2, y=value*100, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
    geom_point(data=diag.data[(3*length.diag+1):(4*length.diag),],aes(x=1:length.diag, y=diag.rate*100), color="red", shape=15, size=2) +
    labs(title="All F: 25-39y", x="Year", y="Reported diagnoses per 100") +
    coord_cartesian(ylim=c(0,max.o.f))+
    scale_x_discrete(labels=seq((end.year-length.diag+1),end.year,1))

  ### plot diagnosed cases by subpop, age, and sex ###
  plot.diag.y.m.1 <- ggplot(data=out.diag.y.m.1)+
    geom_line(aes(x=Var2, y=value*100, group=Var1), color="grey") +
    stat_summary(aes(x=Var2, y=value*100, group=1), geom="line", fun.y="median", color="black", size=1) +
    #geom_line(data=mean.diag.y.m.1, aes(x=1:5, y=mean.diag.y.m.1[,1]), size=1)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=6)) +
    geom_point(data=diag.subpop.rate,aes(x=1:length.rr, y=y.m.1*100), color="red", shape=15, size=2) +
    labs(title="Black M (+MSM): 15-24y", x="Year", y="Reported diagnoses per 100") +
    coord_cartesian(ylim=c(0,max.y.m))+
    scale_x_discrete(labels=seq((end.year-length.rr+1),end.year,1))

  plot.diag.o.m.1 <- ggplot(data=out.diag.o.m.1)+
    geom_line(aes(x=Var2, y=value*100, group=Var1),color="grey") +
    stat_summary(aes(x=Var2, y=value*100, group=1), geom="line", fun.y="median", color="black", size=1) +
    #geom_line(data=mean.diag.o.m.1, aes(x=1:5, y=mean.diag.o.m.1[,1]), size=1)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=6)) +
    geom_point(data=diag.subpop.rate,aes(x=1:length.rr, y=o.m.1*100), color="red", shape=15, size=2) +
    labs(title="Black M (+MSM): 25-39y", x="Year", y="Reported diagnoses per 100") +
    coord_cartesian(ylim=c(0,max.o.m))+
    scale_x_discrete(labels=seq((end.year-length.rr+1),end.year,1))

  plot.diag.y.m.2<- ggplot(data=out.diag.y.m.2)+
    geom_line(aes(x=Var2, y=value*100, group=Var1), color="grey") +
    stat_summary(aes(x=Var2, y=value*100, group=1), geom="line", fun.y="median", color="black", size=1) +
    #geom_line(data=mean.diag.y.m.2, aes(x=1:5, y=mean.diag.y.m.2[,1]), size=1)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=6)) +
    geom_point(data=diag.subpop.rate,aes(x=1:length.rr, y=100*y.m.2), color="red", shape=15, size=2) +
    labs(title="Other M (+MSM): 15-24y", x="Year", y="Reported diagnoses per 100") +
    coord_cartesian(ylim=c(0,max.y.m))+
    scale_x_discrete(labels=seq((end.year-length.rr+1),end.year,1))

  plot.diag.o.m.2<- ggplot(data=out.diag.o.m.2)+
    geom_line(aes(x=Var2, y=value*100, group=Var1),color="grey") +
    stat_summary(aes(x=Var2, y=value*100, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=6)) +
    geom_point(data=diag.subpop.rate,aes(x=1:length.rr, y=100*o.m.2), color="red", shape=15, size=2) +
    labs(title="Other M (+MSM): 25-39y", x="Year",y="Reported diagnoses per 100") +
    coord_cartesian(ylim=c(0,max.o.m))+
    scale_x_discrete(labels=seq((end.year-length.rr+1),end.year,1))

  plot.diag.y.m.3<- ggplot(data=out.diag.y.m.3)+
    geom_line(aes(x=Var2, y=value*100, group=Var1),color="grey") +
    stat_summary(aes(x=Var2, y=value*100, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=6)) +
    geom_point(data=diag.subpop.rate,aes(x=1:length.rr, y=100*y.m.3), color="red", shape=15, size=2) +
    labs(title="Hispanic M (+MSM): 15-24y", x="Year", y="Reported diagnoses per 100") +
    coord_cartesian(ylim=c(0,max.y.m))+
    scale_x_discrete(labels=seq((end.year-length.rr+1),end.year,1))

  plot.diag.o.m.3<- ggplot(data=out.diag.o.m.3)+
    geom_line(aes(x=Var2, y=value*100, group=Var1),color="grey") +
    stat_summary(aes(x=Var2, y=value*100, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=6)) +
    geom_point(data=diag.subpop.rate,aes(x=1:length.rr, y=100*o.m.3), color="red", shape=15, size=2) +
    labs(title="Hispanic M (+MSM): 25-39y", x="Year", y="Reported diagnoses per 100") +
    coord_cartesian(ylim=c(0,max.o.m))+
    scale_x_discrete(labels=seq((end.year-length.rr+1),end.year,1))

  plot.diag.y.f.1<- ggplot(data=out.diag.y.f.1)+
    geom_line(aes(x=Var2, y=value*100, group=Var1),color="grey") +
    stat_summary(aes(x=Var2, y=value*100, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=6)) +
    geom_point(data=diag.subpop.rate,aes(x=1:length.rr, y=100*y.f.1), color="red", shape=15, size=2) +
    labs(title="Black F: 15-24y", x="Year", y="Reported diagnoses per 100") +
    coord_cartesian(ylim=c(0,max.y.f))+
    scale_x_discrete(labels=seq((end.year-length.rr+1),end.year,1))

  plot.diag.o.f.1 <- ggplot(data=out.diag.o.f.1)+
    geom_line(aes(x=Var2, y=value*100, group=Var1),color="grey") +
    stat_summary(aes(x=Var2, y=value*100, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=6)) +
    geom_point(data=diag.subpop.rate,aes(x=1:length.rr, y=100*o.f.1), color="red", shape=15, size=2) +
    labs(title="Black F: 25-39y", x="Year", y="Reported diagnoses per 100") +
    coord_cartesian(ylim=c(0,max.o.f))+
    scale_x_discrete(labels=seq((end.year-length.rr+1),end.year,1))

  plot.diag.y.f.2<- ggplot(data=out.diag.y.f.2)+
    geom_line(aes(x=Var2, y=value*100, group=Var1),color="grey") +
    stat_summary(aes(x=Var2, y=value*100, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=6)) +
    geom_point(data=diag.subpop.rate,aes(x=1:length.rr, y=100*y.f.2), color="red", shape=15, size=2) +
    labs(title="Other F: 15-24y", x="Year", y="Reported diagnoses per 100") +
    coord_cartesian(ylim=c(0,max.y.f))+
    scale_x_discrete(labels=seq((end.year-length.rr+1),end.year,1))

  plot.diag.o.f.2<- ggplot(data=out.diag.o.f.2)+
    geom_line(aes(x=Var2, y=value*100, group=Var1),color="grey") +
    stat_summary(aes(x=Var2, y=value*100, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=6)) +
    geom_point(data=diag.subpop.rate,aes(x=1:length.rr, y=100*o.f.2), color="red", shape=15, size=2) +
    labs(title="Other F: 25-39y", x="Year",y="Reported diagnoses per 100") +
    coord_cartesian(ylim=c(0,max.o.f))+
    scale_x_discrete(labels=seq((end.year-length.rr+1),end.year,1))

  plot.diag.y.f.3 <- ggplot(data=out.diag.y.f.3)+
    geom_line(aes(x=Var2, y=value*100, group=Var1),color="grey") +
    stat_summary(aes(x=Var2, y=value*100, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=6)) +
    geom_point(data=diag.subpop.rate,aes(x=1:length.rr, y=100*y.f.3), color="red", shape=15, size=2) +
    labs(title="Hispanic F: 15-24y", x="Year", y="Reported diagnoses per 100") +
    coord_cartesian(ylim=c(0,max.y.f))+
    scale_x_discrete(labels=seq((end.year-length.rr+1),end.year,1))

  plot.diag.o.f.3 <- ggplot(data=out.diag.o.f.3)+
    geom_line(aes(x=Var2, y=value*100, group=Var1),color="grey") +
    stat_summary(aes(x=Var2, y=value*100, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=6)) +
    geom_point(data=diag.subpop.rate,aes(x=1:length.rr, y=100*o.f.3), color="red", shape=15, size=2) +
    labs(title="Hispanic F: 25-39y", x="Year", y="Reported diagnoses per 100") +
    coord_cartesian(ylim=c(0,max.o.f))+
    scale_x_discrete(labels=seq((end.year-length.rr+1),end.year,1))

  plot.diag.y.msm <- ggplot(data=out.diag.y.msm)+
    geom_line(aes(x=Var2, y=value*100, group=Var1), color="grey") +
    stat_summary(aes(x=Var2, y=value*100, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=6)) +
    #geom_point(data=diag.subpop.data[1:15,],aes(x=1:15, y=diag.sub.time.dat), color="red", shape=15, size=2) +
    labs(title="MSM: 15-24y", x="Year", y="Reported diagnoses per 100") +
    coord_cartesian(ylim=c(0,max.msm))+
    scale_x_discrete(labels=seq((end.year-length.rr+1),end.year,1))

  plot.diag.o.msm <- ggplot(data=out.diag.o.msm)+
    geom_line(aes(x=Var2, y=value*100, group=Var1), color="grey") +
    stat_summary(aes(x=Var2, y=value*100, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=6)) +
    #geom_point(data=diag.subpop.data[1:15,],aes(x=1:15, y=diag.sub.time.dat), color="red", shape=15, size=2) +
    labs(title="MSM: 25-39y", x="Year", y="Reported diagnoses per 100") +
    coord_cartesian(ylim=c(0,max.msm))+
    scale_x_discrete(labels=seq((end.year-length.rr+1),end.year,1))

  #get max values of incidence for setting y-axis in plots
  max.inc.y.m <- max(out.inc.y.m.1$value, out.inc.y.m.2$value, out.inc.y.m.3$value, na.rm=TRUE)
  max.inc.o.m <- max(out.inc.o.m.1$value, out.inc.o.m.2$value, out.inc.o.m.3$value, na.rm=TRUE)
  max.inc.y.f <-max(out.inc.y.f.1$value, out.inc.y.f.2$value, out.inc.y.f.3$value, na.rm=TRUE)
  max.inc.o.f <-max(out.inc.o.f.1$value, out.inc.o.f.2$value, out.inc.o.f.3$value, na.rm=TRUE)
  max.inc <- ceiling(max(max.inc.y.m, max.inc.o.m, max.inc.y.f, max.inc.o.f))
  max.inc.msm <-ceiling(max(out.inc.y.msm$value, out.inc.o.msm$value, na.rm=TRUE)+1)

  # out.inc.msm <- melt(100*as.matrix(subset(pred, select=(prev.inc1+6*length.diag):(prev.inc1+8*length.diag-1)))/sum(n.i[m4]))

  # out.inc.msm <- melt(100*as.matrix(subset(pred, select=(prev.inc1+6*length.diag):(prev.inc1+7*length.diag-1))))
  # out.inc.msm$value <- out.inc.msm$value +
    # melt(100*as.matrix(subset(pred, select=(prev.inc1+7*length.diag):(prev.inc1+8*length.diag-1))))$value
  out.inc.msm <- out.inc.y.msm
  out.inc.msm$value <- (out.inc.y.msm$value * sum(n.i[y.m.4])) + (out.inc.o.msm$value * sum(n.i[o.m.4]))
  out.inc.msm$value <- out.inc.msm$value / sum(n.i[m4])

  out.inc.msm %>% group_by(Var1) %>%
    summarize(mean = mean(value)) -> out.inc.msm.mean

  plot.inc.msm.mean <-
    ggplot(data = out.inc.msm.mean, aes(x = mean, y = ..ncount..)) +
    geom_histogram(
      bins = 30,
      fill = "cornflowerblue", size = 0.1, colour = "cornflowerblue", alpha =
        0.6
    ) +
    geom_area(
      data = data.frame(x = seq(0, 50, 0.1), y = dgamma(seq(0, 50, 0.1), 3.358014, 0.389351)),
      mapping = aes(x = x, y = y / max(y)),
      fill = 'dimgrey', alpha = 0.3) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    ylab("density") +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    expand_limits(x=c(0.5,1))+
    labs(x=paste0("Mean MSM Incidence During 2010-", gc_env$end.year))



  ### plot incident cases by subpop, age, and sex ###

  plot.inc.y.m.1 <- ggplot(data=out.inc.y.m.1)+
    geom_line(aes(x=Var2, y=value, group=Var1), color="grey") +
    stat_summary(aes(x=Var2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
    labs(title="Black M: 15-24y", x="Year", y="Incidence (%)") +
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-length.diag+1),end.year,1))

  plot.inc.o.m.1<- ggplot(data=out.inc.o.m.1)+
    geom_line(aes(x=Var2, y=value, group=Var1),color="grey") +
    stat_summary(aes(x=Var2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
    labs(title="Black M: 25-39y", x="Year", y="Incidence (%)") +
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-length.diag+1),end.year,1))

  plot.inc.y.m.2<- ggplot(data=out.inc.y.m.2)+
    geom_line(aes(x=Var2, y=value, group=Var1), color="grey") +
    stat_summary(aes(x=Var2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
    labs(title="Other M: 15-24y", x="Year", y="Incidence (%)") +
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-length.diag+1),end.year,1))

  plot.inc.o.m.2<- ggplot(data=out.inc.o.m.2)+
    geom_line(aes(x=Var2, y=value, group=Var1),color="grey") +
    stat_summary(aes(x=Var2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
    labs(title="Other M: 25-39y", x="Year", y="Incidence (%)") +
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-length.diag+1),end.year,1))

  plot.inc.y.m.3<- ggplot(data=out.inc.y.m.3)+
    geom_line(aes(x=Var2, y=value, group=Var1),color="grey") +
    stat_summary(aes(x=Var2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
    labs(title="Hispanic M: 15-24y", x="Year", y="Incidence (%)") +
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-length.diag+1),end.year,1))

  plot.inc.o.m.3<- ggplot(data=out.inc.o.m.3)+
    geom_line(aes(x=Var2, y=value, group=Var1),color="grey") +
    stat_summary(aes(x=Var2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
    labs(title="Hispanic M: 25-39y", x="Year", y="Incidence (%)") +
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-length.diag+1),end.year,1))

  plot.inc.y.f.1<- ggplot(data=out.inc.y.f.1)+
    geom_line(aes(x=Var2, y=value, group=Var1),color="grey") +
    stat_summary(aes(x=Var2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
    labs(title="Black F: 15-24y", x="Year", y="Incidence (%)") +
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-length.diag+1),end.year,1))

  plot.inc.o.f.1<- ggplot(data=out.inc.o.f.1)+
    geom_line(aes(x=Var2, y=value, group=Var1),color="grey") +
    stat_summary(aes(x=Var2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
    labs(title="Black F: 25-39y", x="Year", y="Incidence (%)") +
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-length.diag+1),end.year,1))

  plot.inc.y.f.2<- ggplot(data=out.inc.y.f.2)+
    geom_line(aes(x=Var2, y=value, group=Var1),color="grey") +
    stat_summary(aes(x=Var2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
    labs(title="Other F: 15-24y", x="Year", y="Incidence (%)") +
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-length.diag+1),end.year,1))

  plot.inc.o.f.2<- ggplot(data=out.inc.o.f.2)+
    geom_line(aes(x=Var2, y=value, group=Var1),color="grey") +
    stat_summary(aes(x=Var2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
    labs(title="Other F: 25-39y", x="Year", y="Incidence (%)") +
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-length.diag+1),end.year,1))

  plot.inc.y.f.3<- ggplot(data=out.inc.y.f.3)+
    geom_line(aes(x=Var2, y=value, group=Var1),color="grey") +
    stat_summary(aes(x=Var2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
    labs(title="Hispanic F: 15-24y", x="Year", y="Incidence (%)") +
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-length.diag+1),end.year,1))

  plot.inc.o.f.3 <- ggplot(data=out.inc.o.f.3)+
    geom_line(aes(x=Var2, y=value, group=Var1),color="grey") +
    stat_summary(aes(x=Var2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
    labs(title="Hispanic F: 25-39y", x="Year", y="Incidence (%)") +
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-length.diag+1),end.year,1))

  plot.inc.y.msm <- ggplot(data=out.inc.y.msm)+
    geom_line(aes(x=Var2, y=value, group=Var1), color="grey") +
    stat_summary(aes(x=Var2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
    labs(title="MSM: 15-24y", x="Year", y="Incidence (%)") +
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-length.diag+1),end.year,1))

  plot.inc.o.msm <- ggplot(data=out.inc.o.msm)+
    geom_line(aes(x=Var2, y=value, group=Var1), color="grey") +
    stat_summary(aes(x=Var2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
    labs(title="MSM: 25-39y", x="Year", y="Incidence (%)") +
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-length.diag+1),end.year,1))

  #get max values of rr reported cases for setting y-axis
  max.rr.diag <- max(out.rr.diag.m.y.1$value, out.rr.diag.m.o.1$value, out.rr.diag.m.y.2$value, out.rr.diag.m.o.2$value,
                     out.rr.diag.f.y.1$value, out.rr.diag.f.o.1$value, out.rr.diag.f.y.2$value, out.rr.diag.f.o.2$value, diag.rr$rr_diag*1.2, na.rm=TRUE)
  max.rr.diag.1 <- max(out.rr.diag.m.y.1$value, out.rr.diag.m.o.1$value,
                       out.rr.diag.f.y.1$value, out.rr.diag.f.o.1$value, diag.rr$rr_diag*1.2, na.rm=TRUE)+1
  max.rr.diag.2 <- max(out.rr.diag.m.y.2$value, out.rr.diag.m.o.2$value,
                       out.rr.diag.f.y.2$value, out.rr.diag.f.o.2$value, diag.rr$rr_diag*1.2, na.rm=TRUE)+0.5

  plot.rr.diag.m.y.1 <- ggplot(data=out.rr.diag.m.y.1)+
    geom_line(aes(x=Var2, y=value, group=Var1), color="grey") +
    stat_summary(aes(x=Var2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    geom_pointrange(data=diag.rr[1:length.rr,],aes(x=1:length.rr, y=rr_diag, ymin=rr_diag*0.8, ymax=rr_diag*1.2), color="red", shape=15, size=0.5) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
    labs(title="Black vs all M: 15-24", x="Year", y="Reported diagnoses relative risk") +
    coord_cartesian(ylim=c(0,max.rr.diag.1))+
    scale_x_discrete(labels=seq((end.year-length.rr+1),end.year,1))

  plot.rr.diag.m.o.1 <- ggplot(data=out.rr.diag.m.o.1)+
    geom_line(aes(x=Var2, y=value, group=Var1), color="grey") +
    stat_summary(aes(x=Var2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    geom_pointrange(data=diag.rr[(length.rr+1):(2*length.rr),],aes(x=1:length.rr, y=rr_diag, ymin=rr_diag*0.8, ymax=rr_diag*1.2), color="red", shape=15, size=0.5) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
    labs(title="Black vs all M: 25-39", x="Year", y="Reported diagnoses relative risk") +
    coord_cartesian(ylim=c(0,max.rr.diag.1))+
    scale_x_discrete(labels=seq((end.year-length.rr+1),end.year,1))

  plot.rr.diag.m.y.2 <- ggplot(data=out.rr.diag.m.y.2)+
    geom_line(aes(x=Var2, y=value, group=Var1), color="grey") +
    stat_summary(aes(x=Var2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    geom_pointrange(data=diag.rr[(2*length.rr+1):(3*length.rr),],aes(x=1:length.rr, y=rr_diag, ymin=rr_diag*0.8, ymax=rr_diag*1.2), color="red", shape=15, size=0.5) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
    labs(title="Hispanic vs all M: 15-24", x="Year", y="Reported diagnoses relative risk") +
    coord_cartesian(ylim=c(0,max.rr.diag.2))+
    scale_x_discrete(labels=seq((end.year-length.rr+1),end.year,1))

  plot.rr.diag.m.o.2 <- ggplot(data=out.rr.diag.m.o.2)+
    geom_line(aes(x=Var2, y=value, group=Var1), color="grey") +
    stat_summary(aes(x=Var2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    geom_pointrange(data=diag.rr[(3*length.rr+1):(4*length.rr),],aes(x=1:length.rr, y=rr_diag, ymin=rr_diag*0.8, ymax=rr_diag*1.2), color="red", shape=15, size=0.5) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
    labs(title="Hispanic vs all M: 25-39", x="Year", y="Reported diagnoses relative risk") +
    coord_cartesian(ylim=c(0,max.rr.diag.2))+
    scale_x_discrete(labels=seq((end.year-length.rr+1),end.year,1))

  plot.rr.diag.f.y.1 <- ggplot(data=out.rr.diag.f.y.1)+
    geom_line(aes(x=Var2, y=value, group=Var1), color="grey") +
    stat_summary(aes(x=Var2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    geom_pointrange(data=diag.rr[(4*length.rr+1):(5*length.rr),],aes(x=1:length.rr, y=rr_diag, ymin=rr_diag*0.8, ymax=rr_diag*1.2), color="red", shape=15, size=0.5) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
    labs(title="Black vs all F: 15-24", x="Year", y="Reported diagnoses relative risk") +
    coord_cartesian(ylim=c(0,max.rr.diag.1))+
    scale_x_discrete(labels=seq((end.year-length.rr+1),end.year,1))

  plot.rr.diag.f.o.1 <- ggplot(data=out.rr.diag.f.o.1)+
    geom_line(aes(x=Var2, y=value, group=Var1), color="grey") +
    stat_summary(aes(x=Var2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    geom_pointrange(data=diag.rr[(5*length.rr+1):(6*length.rr),],aes(x=1:length.rr, y=rr_diag, ymin=rr_diag*0.8, ymax=rr_diag*1.2), color="red", shape=15, size=0.5) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
    labs(title="Black vs all F: 25-39", x="Year", y="Reported diagnoses relative risk") +
    coord_cartesian(ylim=c(0,max.rr.diag.1))+
    scale_x_discrete(labels=seq((end.year-length.rr+1),end.year,1))

  plot.rr.diag.f.y.2 <- ggplot(data=out.rr.diag.f.y.2)+
    geom_line(aes(x=Var2, y=value, group=Var1), color="grey") +
    stat_summary(aes(x=Var2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    #geom_line(data=mean.rr.y.b, aes(x=1:12, y=mean.rr.y.b[,1]), size=1)+
    geom_pointrange(data=diag.rr[(6*length.rr+1):(7*length.rr),],aes(x=1:length.rr, y=rr_diag, ymin=rr_diag*0.8, ymax=rr_diag*1.2), color="red", shape=15, size=0.5) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
    labs(title="Hispanic vs all F: 15-24", x="Year", y="Reported diagnoses relative risk") +
    coord_cartesian(ylim=c(0,max.rr.diag.2))+
    scale_x_discrete(labels=seq((end.year-length.rr+1),end.year,1))

  plot.rr.diag.f.o.2 <- ggplot(data=out.rr.diag.f.o.2)+
    geom_line(aes(x=Var2, y=value, group=Var1), color="grey") +
    stat_summary(aes(x=Var2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    geom_pointrange(data=diag.rr[(7*length.rr+1):(8*length.rr),],aes(x=1:length.rr, y=rr_diag, ymin=rr_diag*0.8, ymax=rr_diag*1.2), color="red", shape=15, size=0.5) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
    labs(title="Hispanic vs all F: 25-39", x="Year", y="Reported diagnoses relative risk") +
    coord_cartesian(ylim=c(0, max.rr.diag.2))+
    scale_x_discrete(labels=seq((end.year-length.rr+1),end.year,1))

  ### Additional plots ###

  plot.subpop <- ggplot(data=out.subpop)+
    geom_boxplot(aes(x=Var2, y=value, fill=Var2), outlier.alpha = 0.3) +
    geom_pointrange(data=s.dist.sd,aes(x=1:nrow(s.dist.sd), y=dat.s.dist,ymin=lcl,ymax=ucl), color="dimgrey", shape=15, size=0.5, alpha=0.8) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1),title=element_text(size=8) ) +
    labs(title="Subpopulation assortative mixing\n(not used in calibration)", x="Population", y="Proportion of partnerships \nwith same subpopulation") +
    coord_cartesian(ylim=c(0,1))+
    scale_x_discrete(labels=c("M black","M other","M Hispanic", "F black","F other", "F Hispanic"))

  plot.age <- ggplot(data=out.age)+
    geom_boxplot(aes(x=Var2, y=value, fill=Var2), outlier.alpha = 0.3) +
    geom_pointrange(data=age.data,aes(x=1:nrow(age.data), y=p.same.age, ymin=lcl, ymax=ucl), color="red", shape=15, size=0.5) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1),title=element_text(size=8)) +
    labs(title="Age assortative mixing", x="Population", y="Proportion of partnerships \nwith same age group") +
    coord_cartesian(ylim=c(0,1))+
    scale_x_discrete(labels=c("M young","M old","F young", "F old"))

  color_scale_psymp <- if (gc_env$site == 'SF') c('red',rep('dimgrey', length=nrow(out.symp)-1)) else rep('red', length=nrow(out.symp))
  shape_scale_psymp <- if (gc_env$site == 'SF') c(15, rep(4, length = nrow(out.symp)-1)) else rep(15, length=nrow(out.symp))

  plot.p.symp <- ggplot(data=out.symp)+
    geom_boxplot(aes(x=Var2, y=value, fill=Var2), outlier.alpha = 0.3) +
    geom_pointrange(data=p.symp.dat,
                    aes(x=1:length(p.symp.ssun), y=p.symp,
                        ymax = p.symp + sqrt(p.symp.var),
                        ymin = p.symp - sqrt(p.symp.var),
                        color = as.factor(ifelse(sex == 'M' & age_cat==1, 1, 0)),
                        shape = as.factor(ifelse(sex == 'M' & age_cat==1, 1, 0))),
                    size=.5) +
    theme_classic() +
  theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
  labs(title="Proportion of diagnoses\nthat are symptomatic", x="Population", y="Proportion of diagnoses that are symptomatic") +
  coord_cartesian(ylim=c(0,1))+
  scale_x_discrete(labels=c("MSW 15-24 y","MSW 25-39 y", "MSM 15-24 y", "MSM 25-39 y","F 15-24 y", "F 25-39 y")) +
  scale_color_manual(values=color_scale_psymp) +
  scale_shape_manual(values=shape_scale_psymp)


  plot.p.msm.y <- ggplot(data=out.p.msm.y)+
    geom_line(aes(x=Var2, y = value, group = Var1), color = 'grey') +
    stat_summary(fun.y = mean, geom='line', color = 'black', aes(x = Var2, y=value, group=1)) +
    geom_point(data=as.data.frame(p.msm.plot[1:length.msm]),aes(x=1:length.msm, y=p.msm.plot[1:length.msm]), color="red", shape=15, size=2) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
    labs(title="Proportion of male diagnoses\n in MSM 15-24 y", x="Year", y="Proportion of diagnoses in MSM") +
    coord_cartesian(ylim=c(0,1))+
    scale_x_discrete(labels=seq(end.year-length.msm + 1, end.year, 1))

  plot.p.msm.o <- ggplot(data=out.p.msm.o)+
    geom_line(aes(x=Var2, y = value, group = Var1), color = 'grey') +
    stat_summary(fun.y = mean, geom='line', color = 'black', aes(x = Var2, y=value, group=1)) +
    geom_point(data=as.data.frame(p.msm.plot[(length.msm+1):(2*length.msm)]), aes(x=1:length.msm, y=p.msm.plot[(length.msm+1):(2*length.msm)]), color="red", shape=15, size=2) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
    labs(title="Proportion of male diagnoses\n in MSM 25-39y y", x="Year", y="Proportion of diagnoses in MSM") +
    coord_cartesian(ylim=c(0,1))+
    scale_x_discrete(labels=seq((end.year-length.msm+1),end.year,1))

  ######################################
  ###   Plot priors and posteriors   ###
  ######################################

  #epsilon.1
  epsilon.1.post<-data.frame(var1=ilogit(trace.burn.thin[,"logit.epsilon.1"]))
  epsilon.1.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["epsilon.1"],prior.param2["epsilon.1"])))
  plot.epsilon.1 <- ggplot() +
    geom_histogram(data=epsilon.1.post, aes(x=var1, y=..ncount..), fill="slateblue",size=0.1, colour="slateblue", alpha=0.6, binwidth=0.025)+
    geom_area(data=epsilon.1.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    ylab("density") +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    expand_limits(x=c(0,1))+
    labs(x="Risk group mixing:\n black")

  #epsilon.2
  epsilon.2.post<-data.frame(var1=ilogit(trace.burn.thin[,"logit.epsilon.2"]))
  epsilon.2.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["epsilon.2"],prior.param2["epsilon.2"])))
  plot.epsilon.2 <- ggplot() +
    geom_histogram(data=epsilon.2.post, aes(x=var1, y=..ncount..), fill="slateblue",size=0.1, colour="slateblue", alpha=0.6, binwidth=0.025)+
    geom_area(data=epsilon.2.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    ylab("density") +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    expand_limits(x=c(0,1))+
    labs(x="Risk group mixing:\n other")

  #epsilon.3
  epsilon.3.post<-data.frame(var1=ilogit(trace.burn.thin[,"logit.epsilon.3"]))
  epsilon.3.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["epsilon.3"],prior.param2["epsilon.3"])))
  plot.epsilon.3 <- ggplot() +
    geom_histogram(data=epsilon.3.post, aes(x=var1, y=..ncount..), fill="slateblue",size=0.1, colour="slateblue", alpha=0.6, binwidth=0.025)+
    geom_area(data=epsilon.3.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    ylab("density") +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    expand_limits(x=c(0,1))+
    labs(x="Risk group mixing:\n Hispanic")

  #epsilon.4
  epsilon.4.post<-data.frame(var1=ilogit(trace.burn.thin[,"logit.epsilon.4"]))
  epsilon.4.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["epsilon.4"],prior.param2["epsilon.4"])))
  plot.epsilon.4 <- ggplot() +
    geom_histogram(data=epsilon.4.post, aes(x=var1, y=..ncount..), fill="slateblue",size=0.1, colour="slateblue", alpha=0.6, binwidth=0.025)+
    geom_area(data=epsilon.4.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    ylab("density") +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    expand_limits(x=c(0,1))+
    labs(x="Risk group mixing:\n MSM")

  #pi.m
  pi.m.post<-data.frame(var1=ilogit(trace.burn.thin[,"logit.pi.m"]))
  pi.m.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["pi.m"],prior.param2["pi.m"])))
  plot.pi.m <- ggplot() +
    geom_histogram(data=pi.m.post, aes(x=var1, y=..ncount..), fill="dodgerblue",size=0.1, colour="dodgerblue", alpha=0.6, binwidth=0.025)+
    geom_area(data=pi.m.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    ylab("density") +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    expand_limits(x=c(0,1))+
    labs(x="Age assortativity:\n young M/old F")

  #pi.f
  pi.f.post<-data.frame(var1=ilogit(trace.burn.thin[,"logit.pi.f"]))
  pi.f.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["pi.f"],prior.param2["pi.f"])))
  plot.pi.f <- ggplot() +
    geom_histogram(data=pi.f.post, aes(x=var1, y=..ncount..), fill="dodgerblue",size=0.1, colour="dodgerblue", alpha=0.6, binwidth=0.025)+
    geom_area(data=pi.f.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    ylab("density") +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    expand_limits(x=c(0,1))+
    labs(x="Age assortativity:\n young F/old M")

  #pi.msm
  pi.msm.post<-data.frame(var1=ilogit(trace.burn.thin[,"logit.pi.msm"]))
  pi.msm.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["pi.msm"],prior.param2["pi.msm"])))
  plot.pi.msm <- ggplot() +
    geom_histogram(data=pi.msm.post, aes(x=var1, y=..ncount..), fill="dodgerblue",size=0.1, colour="dodgerblue", alpha=0.6, binwidth=0.025)+
    geom_area(data=pi.msm.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    ylab("density") +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    expand_limits(x=c(0,1))+
    labs(x="Age assortativity:\n MSM")

  #theta.1
  theta.1.post<-data.frame(var1=ilogit(trace.burn.thin[,"logit.theta.1"]))
  theta.1.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["theta.1"],prior.param2["theta.1"])))
  plot.theta.1 <- ggplot() +
    geom_histogram(data=theta.1.post, aes(x=var1, y=..ncount..), fill="cornflowerblue",size=0.1, colour="cornflowerblue", alpha=0.6, binwidth=0.01)+
    geom_area(data=theta.1.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    ylab("density") +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    expand_limits(x=c(0.5,1))+
    labs(x="Subpopulation mixing:\n black M")

  #theta.2
  theta.2.post<-data.frame(var1=ilogit(trace.burn.thin[,"logit.theta.2"]))
  theta.2.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["theta.2"],prior.param2["theta.2"])))
  plot.theta.2 <- ggplot() +
    geom_histogram(data=theta.2.post, aes(x=var1, y=..ncount..), fill="cornflowerblue",size=0.1, colour="cornflowerblue", alpha=0.6, binwidth=0.01)+
    geom_area(data=theta.2.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    ylab("density") +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    expand_limits(x=c(0.5,1))+
    labs(x="Subpopulation mixing:\n other M")

  #theta.3
  theta.3.post<-data.frame(var1=ilogit(trace.burn.thin[,"logit.theta.3"]))
  theta.3.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["theta.3"],prior.param2["theta.3"])))
  plot.theta.3 <- ggplot() +
    geom_histogram(data=theta.3.post, aes(x=var1, y=..ncount..), fill="cornflowerblue",size=0.1, colour="cornflowerblue", alpha=0.6, binwidth=0.01)+
    geom_area(data=theta.3.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    ylab("density") +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    expand_limits(x=c(0.5,1))+
    labs(x="Subpopulation mixing:\n Hispanic M")

  #theta.4
  theta.4.post<-data.frame(var1=ilogit(trace.burn.thin[,"logit.theta.4"]))
  theta.4.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["theta.4"],prior.param2["theta.4"])))
  plot.theta.4 <- ggplot() +
    geom_histogram(data=theta.4.post, aes(x=var1, y=..ncount..), fill="cornflowerblue",size=0.1, colour="cornflowerblue", alpha=0.6, binwidth=0.01)+
    geom_area(data=theta.4.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    ylab("density") +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    expand_limits(x=c(0.5,1))+
    labs(x="Subpopulation mixing:\n MSM")

  #theta.5
  theta.5.post<-data.frame(var1=ilogit(trace.burn.thin[,"logit.theta.5"]))
  theta.5.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["theta.5"],prior.param2["theta.5"])))
  plot.theta.5 <- ggplot() +
    geom_histogram(data=theta.5.post, aes(x=var1, y=..ncount..), fill="cornflowerblue",size=0.1, colour="cornflowerblue", alpha=0.6, binwidth=0.01)+
    geom_area(data=theta.5.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    ylab("density") +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    expand_limits(x=c(0.5,1))+
    labs(x="Subpopulation mixing:\n black F")

  #theta.6
  theta.6.post<-data.frame(var1=ilogit(trace.burn.thin[,"logit.theta.6"]))
  theta.6.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["theta.6"],prior.param2["theta.6"])))
  plot.theta.6 <- ggplot() +
    geom_histogram(data=theta.6.post, aes(x=var1, y=..ncount..), fill="cornflowerblue",size=0.1, colour="cornflowerblue", alpha=0.6, binwidth=0.01)+
    geom_area(data=theta.6.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    ylab("density") +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    expand_limits(x=c(0.5,1))+
    labs(x="Subpopulation mixing:\n other F")

  #theta.7
  theta.7.post<-data.frame(var1=ilogit(trace.burn.thin[,"logit.theta.7"]))
  theta.7.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["theta.7"],prior.param2["theta.7"])))
  plot.theta.7 <- ggplot() +
    geom_histogram(data=theta.7.post, aes(x=var1, y=..ncount..), fill="cornflowerblue",size=0.1, colour="cornflowerblue", alpha=0.6, binwidth=0.01)+
    geom_area(data=theta.7.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    ylab("density") +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    expand_limits(x=c(0.5,1))+
    labs(x="Subpopulation mixing:\n Hispanic F")

  #c.min.m1
  c.min.m1.post<-data.frame(var1=exp(trace.burn.thin[,"log.c.min.m.1"]))
  c.min.m1.prior <- as.data.frame(cbind(x=seq(0,5,0.01),y=dgamma(seq(0,5,0.01),prior.param1["c.min.m.1"],prior.param2["c.min.m.1"])))
  plot.c.min.m1 <- ggplot() +
    geom_histogram(data=c.min.m1.post, aes(x=var1, y=..ncount..), fill="lightgoldenrod1",size=0.1, colour="lightgoldenrod3", alpha=0.6, binwidth=0.05)+
    geom_area(data=c.min.m1.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=5))) +
    ylab("density") +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    expand_limits(x=c(0,2))+
    labs(x=expression(paste("c"[min], " M, 15-24y")))

  #c.min.m2
  c.min.m2.post<-data.frame(var1=exp(trace.burn.thin[,"log.c.min.m.2"]))
  c.min.m2.prior <- as.data.frame(cbind(x=seq(0,5,0.01),y=dgamma(seq(0,5,0.01),prior.param1["c.min.m.2"],prior.param2["c.min.m.2"])))
  plot.c.min.m2 <- ggplot() +
    geom_histogram(data=c.min.m2.post, aes(x=var1, y=..ncount..), fill="lightgoldenrod1",size=0.1, colour="lightgoldenrod3", alpha=0.6, binwidth=0.05)+
    geom_area(data=c.min.m2.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=5))) +
    ylab("density") +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    expand_limits(x=c(0,2))+
    labs(x=expression(paste("c"[min], " M, 25-39y")))

  #c.min.f1
  c.min.f1.post<-data.frame(var1=exp(trace.burn.thin[,"log.c.min.f.1"]))
  c.min.f1.prior <- as.data.frame(cbind(x=seq(0,5,0.01),y=dgamma(seq(0,5,0.01),prior.param1["c.min.f.1"],prior.param2["c.min.f.1"])))
  plot.c.min.f1 <- ggplot() +
    geom_histogram(data=c.min.f1.post, aes(x=var1, y=..ncount..), fill="lightgoldenrod1",size=0.1, colour="lightgoldenrod3", alpha=0.6, binwidth=0.05)+
    geom_area(data=c.min.f1.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=5))) +
    ylab("density") +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    expand_limits(x=c(0,2))+
    labs(x=expression(paste("c"[min], " F, 15-24y")))

  #c.min.f2
  c.min.f2.post<-data.frame(var1=exp(trace.burn.thin[,"log.c.min.f.2"]))
  c.min.f2.prior <- as.data.frame(cbind(x=seq(0,5,0.01),y=dgamma(seq(0,5,0.01),prior.param1["c.min.f.2"],prior.param2["c.min.f.2"])))
  plot.c.min.f2 <- ggplot() +
    geom_histogram(data=c.min.f2.post, aes(x=var1, y=..ncount..), fill="lightgoldenrod1",size=0.1, colour="lightgoldenrod3", alpha=0.6, binwidth=0.05)+
    geom_area(data=c.min.f2.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=5))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,2))+
    labs(x=expression(paste("c"[min], " F, 25-39y")))

  #c.min.msm1
  c.min.msm1.post<-data.frame(var1=exp(trace.burn.thin[,"log.c.min.msm.1"]))
  c.min.msm1.prior <- as.data.frame(cbind(x=seq(0,5,0.01),y=dgamma(seq(0,5,0.01),prior.param1["c.min.msm.1"],prior.param2["c.min.msm.1"])))
  plot.c.min.msm1 <- ggplot() +
    geom_histogram(data=c.min.msm1.post, aes(x=var1, y=..ncount..), fill="lightgoldenrod1",size=0.1, colour="lightgoldenrod3", alpha=0.6, binwidth=0.05)+
    geom_area(data=c.min.msm1.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=5))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,2))+
    labs(x=expression(paste("c"[min], " MSM, 15-24y")))

  #c.min.msm2
  c.min.msm2.post<-data.frame(var1=exp(trace.burn.thin[,"log.c.min.msm.2"]))
  c.min.msm2.prior <- as.data.frame(cbind(x=seq(0,5,0.01),y=dgamma(seq(0,5,0.01),prior.param1["c.min.msm.2"],prior.param2["c.min.msm.2"])))
  plot.c.min.msm2 <- ggplot() +
    geom_histogram(data=c.min.msm2.post, aes(x=var1, y=..ncount..), fill="lightgoldenrod1",size=0.1, colour="lightgoldenrod3", alpha=0.6, binwidth=0.05)+
    geom_area(data=c.min.msm2.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=5))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,2)) +
    labs(x=expression(paste("c"[min], " MSM, 25-39y")))

  #dur.inf.symp.m
  dur.inf.symp.m.post<-data.frame(var1=exp(trace.burn.thin[,"log.dur.inf.symp.m"]))
  dur.inf.symp.m.prior <- as.data.frame(cbind(x=seq(0,30,0.1),y=dgamma(seq(0,30,0.1),prior.param1["dur.inf.symp.m"],prior.param2["dur.inf.symp.m"])))
  plot.dur.inf.symp.m <- ggplot() +
    geom_histogram(data=dur.inf.symp.m.post, aes(x=var1, y=..ncount..), bins = 30, fill="turquoise",size=0.5, colour="turquoise", alpha=0.6)+
    geom_area(data=dur.inf.symp.m.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,30)) +
    labs(x=expression(paste("Duration symptomatic\ninfection (d): M")))

  #dur.inf.symp.msm
  dur.inf.symp.msm.post<-data.frame(var1=exp(trace.burn.thin[,"log.dur.inf.symp.msm"]))
  dur.inf.symp.msm.prior <- as.data.frame(cbind(x=seq(0,30,0.1),y=dgamma(seq(0,30,0.1),prior.param1["dur.inf.symp.msm"],prior.param2["dur.inf.symp.msm"])))
  plot.dur.inf.symp.msm <- ggplot() +
    geom_histogram(data=dur.inf.symp.msm.post, aes(x=var1, y=..ncount..), bins = 30, fill="turquoise",size=0.5, colour="turquoise", alpha=0.6)+
    geom_area(data=dur.inf.symp.msm.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,30)) +
    labs(x=expression(paste("Duration symptomatic\ninfection (d): MSM")))

  #dur.inf.symp.f
  dur.inf.symp.f.post<-data.frame(var1=exp(trace.burn.thin[,"log.dur.inf.symp.f"]))
  dur.inf.symp.f.prior <- as.data.frame(cbind(x=seq(0,30,0.1),y=dgamma(seq(0,30,0.1),prior.param1["dur.inf.symp.f"],prior.param2["dur.inf.symp.f"])))
  plot.dur.inf.symp.f <- ggplot() +
    geom_histogram(data=dur.inf.symp.f.post, aes(x=var1, y=..ncount..), bins = 30, fill="turquoise",size=0.5, colour="turquoise", alpha=0.6)+
    geom_area(data=dur.inf.symp.f.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    ylab("density") +
    expand_limits(x=c(0, 30)) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    labs(x=expression(paste("Duration symptomatic\ninfection (d): F")))

  #dur.inf.asymp.m
  dur.inf.asymp.m.post<-data.frame(var1=exp(trace.burn.thin[,"log.dur.inf.asymp.m"]))
  dur.inf.asymp.m.prior <- as.data.frame(cbind(x=seq(0,600,1),y=dnorm(seq(0,600,1),prior.param1["dur.inf.asymp.m"],prior.param2["dur.inf.asymp.m"])))
  plot.dur.inf.asymp.m <- ggplot() +
    geom_histogram(data=dur.inf.asymp.m.post, aes(x=var1, y=..ncount..), bins = 30, fill="dodgerblue3",size=0.5, colour="dodgerblue3", alpha=0.6)+
    geom_area(data=dur.inf.asymp.m.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    ylab("density") +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    expand_limits(x=c(0, 600)) +
    labs(x=expression(paste("Duration asymptomatic\ninfection (d): M")))

  #dur.inf.asymp.f
  dur.inf.asymp.f.post<-data.frame(var1=exp(trace.burn.thin[,"log.dur.inf.asymp.f"]))
  dur.inf.asymp.f.prior <- as.data.frame(cbind(x=seq(0,600,1),y=dnorm(seq(0,600,1),prior.param1["dur.inf.asymp.f"],prior.param2["dur.inf.asymp.f"])))
  plot.dur.inf.asymp.f <- ggplot() +
    geom_histogram(data=dur.inf.asymp.f.post, aes(x=var1, y=..ncount..), bins = 30, fill="dodgerblue3",size=0.5, colour="dodgerblue3", alpha=0.6)+
    geom_area(data=dur.inf.asymp.f.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0, 600)) +
    labs(x=expression(paste("Duration asymptomatic\ninfection (d): F")))

  #dur.inf.asymp.msm
  dur.inf.asymp.msm.post<-data.frame(var1=exp(trace.burn.thin[,"log.dur.inf.asymp.msm"]))
  dur.inf.asymp.msm.prior <- as.data.frame(cbind(x=seq(0,600,1),y=dnorm(seq(0,600,1),prior.param1["dur.inf.asymp.msm"],prior.param2["dur.inf.asymp.msm"])))
  plot.dur.inf.asymp.msm <- ggplot() +
    geom_histogram(data=dur.inf.asymp.msm.post, aes(x=var1, y=..ncount..), fill="dodgerblue3",size=0.5, colour="dodgerblue3", alpha=0.6, binwidth=10)+
    geom_area(data=dur.inf.asymp.msm.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,600))+
    labs(x=expression(paste("Duration asymptomatic\ninfection (d): MSM")))

  #b.m
  b.m.post<-data.frame(var1=ilogit(trace.burn.thin[,"logit.b.m"]))
  b.m.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["b.m"],prior.param2["b.m"])))
  plot.b.m <- ggplot() +
    geom_histogram(data=b.m.post, aes(x=var1, y=..ncount..), bins = 30, fill="dodgerblue3",size=0.5, colour="dodgerblue3", alpha=0.6)+
    geom_area(data=b.m.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,1))+
    labs(x=expression(paste("Transmission probability in 2002:\n F to M")))

  #b.f
  b.f.post<-data.frame(var1=ilogit(trace.burn.thin[,"logit.b.f"]))
  b.f.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["b.f"],prior.param2["b.f"])))
  plot.b.f <- ggplot() +
    geom_histogram(data=b.f.post, aes(x=var1, y=..ncount..), fill="dodgerblue3",size=0.5, colour="dodgerblue3", alpha=0.6, binwidth=0.025)+
    geom_area(data=b.f.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    theme(axis.title.x=element_text(size=8))+
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,1))+
    labs(x=expression(paste("Transmission probability in 2002:\nM to F")))

  #b.msm
  b.msm.post<-data.frame(var1=ilogit(trace.burn.thin[,"logit.b.msm"]))
  b.msm.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["b.msm"],prior.param2["b.msm"])))
  plot.b.msm <- ggplot() +
    geom_histogram(data=b.msm.post, aes(x=var1, y=..ncount..), fill="dodgerblue3",size=0.5, colour="dodgerblue3", alpha=0.6, binwidth=0.025)+
    geom_area(data=b.msm.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,1))+
    labs(x=expression(paste("Transmission probability in 2002:\nM to M")))

  #rr.screen.m3
  rr.screen.m3.post<-data.frame(var1=exp(trace.burn.thin[,"log.rr.screen.m3"]))
  rr.screen.m3.prior <- as.data.frame(cbind(x=seq(0,5,0.01),y=dgamma(seq(0,5,0.01),prior.param1["rr.screen.m3"],prior.param2["rr.screen.m3"])))
  plot.rr.screen.m3 <- ggplot() +
    geom_histogram(data=rr.screen.m3.post, aes(x=var1, y=..ncount..), fill="deepskyblue",size=0.5, colour="deepskyblue", alpha=0.6, binwidth=0.1)+
    geom_area(data=rr.screen.m3.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,5))+
    labs(x=expression(paste("RR screen:\nHispanic M")))

  ##rr.screen.f3
  rr.screen.f3.post<-data.frame(var1=exp(trace.burn.thin[,"log.rr.screen.f3"]))
  rr.screen.f3.prior <- as.data.frame(cbind(x=seq(0,5,0.01),y=dgamma(seq(0,5,0.01),prior.param1["rr.screen.f3"],prior.param2["rr.screen.f3"])))
  plot.rr.screen.f3 <- ggplot() +
    geom_histogram(data=rr.screen.f3.post, aes(x=var1, y=..ncount..), fill="deepskyblue",size=0.5, colour="deepskyblue", alpha=0.6, binwidth=0.1)+
    geom_area(data=rr.screen.f3.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,5))+
    labs(x=expression(paste("RR screen:\nHispanic F")))

  ##rr.screen.ac
  rr.screen.ac.post<-data.frame(var1=exp(trace.burn.thin[,"log.rr.screen.ac"]))
  rr.screen.ac.prior <- as.data.frame(cbind(x=seq(0,5,0.01),y=dgamma(seq(0,5,0.01),prior.param1["rr.screen.ac"],prior.param2["rr.screen.ac"])))
  plot.rr.screen.ac <- ggplot() +
    geom_histogram(data=rr.screen.ac.post, aes(x=var1, y=..ncount..), fill="deepskyblue",size=0.5, colour="deepskyblue", alpha=0.6, binwidth=0.1)+
    geom_area(data=rr.screen.ac.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,5))+
    labs(x=expression(paste("RR screen:\nhigh sexual activity group")))

  #symp.m
  symp.m.post<-data.frame(var1=ilogit(trace.burn.thin[,"logit.symp.m"]))
  symp.m.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["symp.m"],prior.param2["symp.m"])))
  plot.symp.m <- ggplot() +
    geom_histogram(data=symp.m.post, aes(x=var1, y=..ncount..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.025)+
    geom_area(data=symp.m.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=5))) +
    ylab("density") +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    expand_limits(x=c(0,1))+
    labs(x="Probability symptomatic:\nM")

  #symp.f
  symp.f.post<-data.frame(var1=ilogit(trace.burn.thin[,"logit.symp.f"]))
  symp.f.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["symp.f"],prior.param2["symp.f"])))
  plot.symp.f <- ggplot() +
    geom_histogram(data=symp.f.post, aes(x=var1, y=..ncount..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.025)+
    geom_area(data=symp.f.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=5))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,1))+
    labs(x="Probability symptomatic:\nF")

  #symp.msm
  symp.msm.post<-data.frame(var1=ilogit(trace.burn.thin[,"logit.symp.msm"]))
  symp.msm.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["symp.msm"],prior.param2["symp.msm"])))
  plot.symp.msm <- ggplot() +
    geom_histogram(data=symp.msm.post, aes(x=var1, y=..ncount..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.025)+
    geom_area(data=symp.msm.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=5))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,1))+
    labs(x="Probability symptomatic:\n MSM")

  #rp.1.1.1.1
  rp.1.1.1.1.post<-data.frame(var1=exp(trace.burn.thin[,"log.rp.1.1.1.1"]))
  rp.1.1.1.1.prior <- as.data.frame(cbind(x=seq(0,10,0.1),y=dgamma(seq(0,10,0.1),prior.param1["rp.1.1.1.1"],prior.param2["rp.1.1.1.1"])))
  plot.rp.1.1.1.1 <- ggplot() +
    geom_histogram(data=rp.1.1.1.1.post, aes(x=var1, y=..ncount..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.5)+
    geom_area(data=rp.1.1.1.1.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,10))+
    labs(x=expression(paste("RR partner change:\nBlack M low AC, 15-24y")))

  #rp.1.1.2.1
  rp.1.1.2.1.post<-data.frame(var1=exp(trace.burn.thin[,"log.rp.1.1.2.1"]))
  rp.1.1.2.1.prior <- as.data.frame(cbind(x=seq(0,100,0.1),y=dnorm(seq(0,100,0.1),prior.param1["rp.1.1.2.1"],prior.param2["rp.1.1.2.1"])))
  plot.rp.1.1.2.1 <- ggplot() +
    geom_histogram(data=rp.1.1.2.1.post, aes(x=var1, y=..ncount..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=5)+
    geom_area(data=rp.1.1.2.1.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,100))+
    labs(x=expression(paste("RR partner change:\nBlack M high AC, 15-24y")))

  #rp.1.2.1.1
  rp.1.2.1.1.post<-data.frame(var1=exp(trace.burn.thin[,"log.rp.1.2.1.1"]))
  rp.1.2.1.1.prior <- as.data.frame(cbind(x=seq(0,10,0.1),y=dgamma(seq(0,10,0.1),prior.param1["rp.1.2.1.1"],prior.param2["rp.1.2.1.1"])))
  plot.rp.1.2.1.1 <- ggplot() +
    geom_histogram(data=rp.1.2.1.1.post, aes(x=var1, y=..ncount..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.5)+
    geom_area(data=rp.1.2.1.1.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    ylab("density") +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    expand_limits(x=c(0,10))+
    labs(x=expression(paste("RR partner change:\nBlack F low AC, 15-24y")))

  #rp.1.2.2.1
  rp.1.2.2.1.post<-data.frame(var1=exp(trace.burn.thin[,"log.rp.1.2.2.1"]))
  rp.1.2.2.1.prior <- as.data.frame(cbind(x=seq(0,40,0.1),y=dgamma(seq(0,40,0.1),prior.param1["rp.1.2.2.1"],prior.param2["rp.1.2.2.1"])))
  plot.rp.1.2.2.1 <- ggplot() +
    geom_histogram(data=rp.1.2.2.1.post, aes(x=var1, y=..ncount..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=2)+
    geom_area(data=rp.1.2.2.1.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,40))+
    labs(x=expression(paste("RR partner change:\nBlack F high AC, 15-24y")))

  #rp.1.1.1.2
  rp.1.1.1.2.post<-data.frame(var1=exp(trace.burn.thin[,"log.rp.1.1.1.2"]))
  rp.1.1.1.2.prior <- as.data.frame(cbind(x=seq(0,5,0.1),y=dgamma(seq(0,5,0.1),prior.param1["rp.1.1.1.2"],prior.param2["rp.1.1.1.2"])))
  plot.rp.1.1.1.2 <- ggplot() +
    geom_histogram(data=rp.1.1.1.2.post, aes(x=var1, y=..ncount..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.2)+
    geom_area(data=rp.1.1.1.2.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,5))+
    labs(x=expression(paste("RR partner change:\nBlack M low AC, 25-39y")))

  #rp.1.1.2.2
  rp.1.1.2.2.post<-data.frame(var1=exp(trace.burn.thin[,"log.rp.1.1.2.2"]))
  rp.1.1.2.2.prior <- as.data.frame(cbind(x=seq(0,100,0.1),y=dnorm(seq(0,100,0.1),prior.param1["rp.1.1.2.2"],prior.param2["rp.1.1.2.2"])))
  plot.rp.1.1.2.2 <- ggplot() +
    geom_histogram(data=rp.1.1.2.2.post, aes(x=var1, y=..ncount..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=5)+
    geom_area(data=rp.1.1.2.2.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,100))+
    labs(x=expression(paste("RR partner change:\nBlack M high AC, 25-39y")))

  #rp.1.2.1.2
  rp.1.2.1.2.post<-data.frame(var1=exp(trace.burn.thin[,"log.rp.1.2.1.2"]))
  rp.1.2.1.2.prior <- as.data.frame(cbind(x=seq(0,5,0.1),y=dgamma(seq(0,5,0.1),prior.param1["rp.1.2.1.2"],prior.param2["rp.1.2.1.2"])))
  plot.rp.1.2.1.2 <- ggplot() +
    geom_histogram(data=rp.1.2.1.2.post, aes(x=var1, y=..ncount..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.2)+
    geom_area(data=rp.1.2.1.2.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,5))+
    labs(x=expression(paste("RR partner change:\nBlack F low AC, 25-39y")))

  #rp.1.2.2.2
  rp.1.2.2.2.post<-data.frame(var1=exp(trace.burn.thin[,"log.rp.1.2.2.2"]))
  rp.1.2.2.2.prior <- as.data.frame(cbind(x=seq(0,40,0.1),y=dgamma(seq(0,40,0.1),prior.param1["rp.1.2.2.2"],prior.param2["rp.1.2.2.2"])))
  plot.rp.1.2.2.2 <- ggplot() +
    geom_histogram(data=rp.1.2.2.2.post, aes(x=var1, y=..ncount..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=2)+
    geom_area(data=rp.1.2.2.2.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,40))+
    labs(x=expression(paste("RR partner change:\nBlack F high AC, 25-39y")))

  #rp.2.1.1.1
  rp.2.1.1.1.post<-data.frame(var1=exp(trace.burn.thin[,"log.rp.2.1.1.1"]))
  rp.2.1.1.1.prior <- as.data.frame(cbind(x=seq(0,10,0.1),y=dgamma(seq(0,10,0.1),prior.param1["rp.2.1.1.1"],prior.param2["rp.2.1.1.1"])))
  plot.rp.2.1.1.1 <- ggplot() +
    geom_histogram(data=rp.2.1.1.1.post, aes(x=var1, y=..ncount..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.2)+
    geom_area(data=rp.2.1.1.1.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,10))+
    labs(x=expression(paste("RR partner change:\nOther M low AC, 15-24y")))

  #rp.2.1.2.1
  rp.2.1.2.1.post<-data.frame(var1=exp(trace.burn.thin[,"log.rp.2.1.2.1"]))
  rp.2.1.2.1.prior <- as.data.frame(cbind(x=seq(0,100,0.1),y=dnorm(seq(0,100,0.1),prior.param1["rp.2.1.2.1"],prior.param2["rp.2.1.2.1"])))
  plot.rp.2.1.2.1 <- ggplot() +
    geom_histogram(data=rp.2.1.2.1.post, aes(x=var1, y=..ncount..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=5)+
    geom_area(data=rp.2.1.2.1.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    ylab("density") +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    expand_limits(x=c(0,100))+
    labs(x=expression(paste("RR partner change:\nOther M high AC, 15-24y")))

  #rp.2.2.1.1
  rp.2.2.1.1.post<-data.frame(var1=exp(trace.burn.thin[,"log.rp.2.2.1.1"]))
  rp.2.2.1.1.prior <- as.data.frame(cbind(x=seq(0,10,0.1),y=dgamma(seq(0,10,0.1),prior.param1["rp.2.2.1.1"],prior.param2["rp.2.2.1.1"])))
  plot.rp.2.2.1.1 <- ggplot() +
    geom_histogram(data=rp.1.2.1.1.post, aes(x=var1, y=..ncount..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.5)+
    geom_area(data=rp.2.2.1.1.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,10))+
    labs(x=expression(paste("RR partner change:\nOther F low AC, 15-24y")))

  #rp.2.2.2.1
  rp.2.2.2.1.post<-data.frame(var1=exp(trace.burn.thin[,"log.rp.2.2.2.1"]))
  rp.2.2.2.1.prior <- as.data.frame(cbind(x=seq(0,40,0.1),y=dgamma(seq(0,40,0.1),prior.param1["rp.2.2.2.1"],prior.param2["rp.2.2.2.1"])))
  plot.rp.2.2.2.1 <- ggplot() +
    geom_histogram(data=rp.2.2.2.1.post, aes(x=var1, y=..ncount..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=1)+
    geom_area(data=rp.2.2.2.1.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,40))+
    labs(x=expression(paste("RR partner change:\nOther F high AC, 15-24y")))

  #rp.2.1.1.2
  rp.2.1.1.2.post<-data.frame(var1=exp(trace.burn.thin[,"log.rp.2.1.1.2"]))
  rp.2.1.1.2.prior <- as.data.frame(cbind(x=seq(0,5,0.1),y=dgamma(seq(0,5,0.1),prior.param1["rp.2.1.1.2"],prior.param2["rp.2.1.1.2"])))
  plot.rp.2.1.1.2 <- ggplot() +
    geom_histogram(data=rp.2.1.1.2.post, aes(x=var1, y=..ncount..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.2)+
    geom_area(data=rp.2.1.1.2.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,5))+
    labs(x=expression(paste("RR partner change:\nOther M low AC, 25-39y")))

  #rp.2.1.2.2
  rp.2.1.2.2.post<-data.frame(var1=exp(trace.burn.thin[,"log.rp.1.1.2.2"]))
  rp.2.1.2.2.prior <- as.data.frame(cbind(x=seq(0,100,0.1),y=dnorm(seq(0,100,0.1),prior.param1["rp.2.1.2.2"],prior.param2["rp.2.1.2.2"])))
  plot.rp.2.1.2.2 <- ggplot() +
    geom_histogram(data=rp.2.1.2.2.post, aes(x=var1, y=..ncount..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=5)+
    geom_area(data=rp.2.1.2.2.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,100))+
    labs(x=expression(paste("RR partner change:\nOther M high AC, 25-39y")))

  #rp.2.2.1.2
  rp.2.2.1.2.post<-data.frame(var1=exp(trace.burn.thin[,"log.rp.2.2.1.2"]))
  rp.2.2.1.2.prior <- as.data.frame(cbind(x=seq(0,5,0.1),y=dgamma(seq(0,5,0.1),prior.param1["rp.2.2.1.2"],prior.param2["rp.2.2.1.2"])))
  plot.rp.2.2.1.2 <- ggplot() +
    geom_histogram(data=rp.2.2.1.2.post, aes(x=var1, y=..ncount..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.2)+
    geom_area(data=rp.2.2.1.2.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,5))+
    labs(x=expression(paste("RR partner change:\nOther F low AC, 25-39y")))

  #rp.2.2.2.2
  rp.2.2.2.2.post<-data.frame(var1=exp(trace.burn.thin[,"log.rp.2.2.2.2"]))
  rp.2.2.2.2.prior <- as.data.frame(cbind(x=seq(0,40,0.1),y=dgamma(seq(0,40,0.1),prior.param1["rp.2.2.2.2"],prior.param2["rp.2.2.2.2"])))
  plot.rp.2.2.2.2 <- ggplot() +
    geom_histogram(data=rp.2.2.2.2.post, aes(x=var1, y=..ncount..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=1)+
    geom_area(data=rp.2.2.2.2.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,40))+
    labs(x=expression(paste("RR partner change:\nOther F high AC, 25-39y")))

  #rp.3.1.2.1
  rp.3.1.2.1.post<-data.frame(var1=exp(trace.burn.thin[,"log.rp.3.1.2.1"]))
  rp.3.1.2.1.prior <- as.data.frame(cbind(x=seq(0,100,0.1),y=dgamma(seq(0,100,0.1),prior.param1["rp.3.1.2.1"],prior.param2["rp.3.1.2.1"])))
  plot.rp.3.1.2.1 <- ggplot() +
    geom_histogram(data=rp.3.1.2.1.post, aes(x=var1, y=..ncount..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=2)+
    geom_area(data=rp.3.1.2.1.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,100))+
    labs(x=expression(paste("RR partner change:\nHispanic M high AC, 15-24y")))

  #rp.3.2.2.1
  rp.3.2.2.1.post<-data.frame(var1=exp(trace.burn.thin[,"log.rp.3.2.2.1"]))
  rp.3.2.2.1.prior <- as.data.frame(cbind(x=seq(0,40,0.1),y=dgamma(seq(0,40,0.1),prior.param1["rp.3.2.2.1"],prior.param2["rp.3.2.2.1"])))
  plot.rp.3.2.2.1 <- ggplot() +
    geom_histogram(data=rp.3.2.2.1.post, aes(x=var1, y=..ncount..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=1)+
    geom_area(data=rp.3.2.2.1.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,40))+
    labs(x=expression(paste("RR partner change:\nHispanic F high AC, 15-24y")))

  #rp.3.1.2.2
  rp.3.1.2.2.post<-data.frame(var1=exp(trace.burn.thin[,"log.rp.3.1.2.2"]))
  rp.3.1.2.2.prior <- as.data.frame(cbind(x=seq(0,100,0.1),y=dgamma(seq(0,100,0.1),prior.param1["rp.3.1.2.2"],prior.param2["rp.3.1.2.2"])))
  plot.rp.3.1.2.2 <- ggplot() +
    geom_histogram(data=rp.3.1.2.2.post, aes(x=var1, y=..ncount..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=2)+
    geom_area(data=rp.3.1.2.2.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,100))+
    labs(x=expression(paste("RR partner change:\nHispanic M high AC, 25-39y")))

  #rp.3.2.2.2
  rp.3.2.2.2.post<-data.frame(var1=exp(trace.burn.thin[,"log.rp.3.2.2.2"]))
  rp.3.2.2.2.prior <- as.data.frame(cbind(x=seq(0,40,0.1),y=dgamma(seq(0,40,0.1),prior.param1["rp.3.2.2.2"],prior.param2["rp.3.2.2.2"])))
  plot.rp.3.2.2.2 <- ggplot() +
    geom_histogram(data=rp.3.2.2.2.post, aes(x=var1, y=..ncount..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=1)+
    geom_area(data=rp.3.2.2.2.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,40))+
    labs(x=expression(paste("RR partner change:\nHispanic F high AC, 25-39y")))

  #rp.4.1.2.1
  rp.4.1.2.1.post<-data.frame(var1=exp(trace.burn.thin[,"log.rp.4.1.2.1"]))
  rp.4.1.2.1.prior <- as.data.frame(cbind(x=seq(0,100,0.1),y=dnorm(seq(0,100,0.1),prior.param1["rp.4.1.2.1"],prior.param2["rp.4.1.2.1"])))
  plot.rp.4.1.2.1 <- ggplot() +
    geom_histogram(data=rp.4.1.2.1.post, aes(x=var1, y=..ncount..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=5)+
    geom_area(data=rp.4.1.2.1.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,100))+
    labs(x=expression(paste("RR partner change:\nMSM high AC, 15-24y")))

  #rp.4.1.2.2
  rp.4.1.2.2.post<-data.frame(var1=exp(trace.burn.thin[,"log.rp.4.1.2.2"]))
  rp.4.1.2.2.prior <- as.data.frame(cbind(x=seq(0,100,0.1),y=dnorm(seq(0,100,0.1),prior.param1["rp.4.1.2.2"],prior.param2["rp.4.1.2.2"])))
  plot.rp.4.1.2.2 <- ggplot() +
    geom_histogram(data=rp.4.1.2.2.post, aes(x=var1, y=..ncount..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=5)+
    geom_area(data=rp.4.1.2.2.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,100))+
    labs(x=expression(paste("RR partner change:\nMSM high AC, 25-39y")))

  #rr.rep.symp.m1
  risk.rep.symp.m1.post<-data.frame(var1=ilogit(trace.burn.thin[,"logit.risk.rep.symp.m1"]))
  risk.rep.symp.m1.prior <- as.data.frame(cbind(x=seq(0,1,0.001),y=dbeta(seq(0,1,0.001),prior.param1["risk.rep.symp.m1"],prior.param2["risk.rep.symp.m1"])))
  plot.rr.rep.symp.m1 <- ggplot() +
    geom_histogram(data=risk.rep.symp.m1.post, aes(x=var1, y=..ncount..), fill="thistle3",size=0.5, colour="thistle3", alpha=0.6, binwidth=0.05)+
    geom_area(data=risk.rep.symp.m1.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8),axis.title.x=element_text(margin=margin(t=10))) +
    ylab("density") +
    expand_limits(x=c(0,1.1))+
    labs(x=expression(paste("RR report symptomatic: \nblack M")))

  #rr.rep.symp.m
  risk.rep.symp.m.post<-data.frame(var1=ilogit(trace.burn.thin[,"logit.risk.rep.symp.m"]))
  risk.rep.symp.m.prior <- as.data.frame(cbind(x=seq(0,1,0.001),y=dbeta(seq(0,1,0.001),prior.param1["risk.rep.symp.m"],prior.param2["risk.rep.symp.m"])))
  plot.rr.rep.symp.m <- ggplot() +
    geom_histogram(data=risk.rep.symp.m.post, aes(x=var1, y=..ncount..), fill="thistle3",size=0.5, colour="thistle3", alpha=0.6, binwidth=0.05)+
    geom_area(data=risk.rep.symp.m.prior,aes(x=x, y=y/max(y)),fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8),axis.title.x=element_text(margin=margin(t=10))) +
    ylab("density") +
    expand_limits(x=c(0,1.1))+
    labs(x=expression(paste("RR report symptomatic:\nHispanic + other M")))

  #rr.rep.symp.f1
  risk.rep.symp.f1.post<-data.frame(var1=ilogit(trace.burn.thin[,"logit.risk.rep.symp.f1"]))
  risk.rep.symp.f1.prior <- as.data.frame(cbind(x=seq(0,1,0.001),y=dbeta(seq(0,1,0.001),prior.param1["risk.rep.symp.f1"],prior.param2["risk.rep.symp.f1"])))
  plot.rr.rep.symp.f1 <- ggplot() +
    geom_histogram(data=risk.rep.symp.f1.post, aes(x=var1, y=..ncount..), fill="thistle3",size=0.5, colour="thistle3", alpha=0.6, binwidth=0.05)+
    geom_area(data=risk.rep.symp.f1.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8),axis.title.x=element_text(margin=margin(t=10))) +
    ylab("density") +
    expand_limits(x=c(0,1.1))+
    labs(x=expression(paste("RR report symptomatic:\nblack F")))

  #rr.rep.symp.f
  risk.rep.symp.f.post<-data.frame(var1=ilogit(trace.burn.thin[,"logit.risk.rep.symp.f"]))
  risk.rep.symp.f.prior <- as.data.frame(cbind(x=seq(0,1,0.001),y=dbeta(seq(0,1,0.001),prior.param1["risk.rep.symp.f"],prior.param2["risk.rep.symp.f"])))
  plot.rr.rep.symp.f <- ggplot() +
    geom_histogram(data=risk.rep.symp.f.post, aes(x=var1, y=..ncount..), fill="thistle3",size=0.5, colour="thistle3", alpha=0.6, binwidth=0.05)+
    geom_area(data=risk.rep.symp.f.prior,aes(x=x, y=y/max(y)), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    ylab("density") +
    expand_limits(x=c(0,1.1))+
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8),axis.title.x=element_text(margin=margin(t=10))) +
    labs(x=expression(paste("RR report symptomatic:\nHispanic + other F")))

  ### prepare time-varying parameters ###

  x<-apply(post.sample$theta,1,update_ctrl)
  post.screen.bc<-sapply(x,`[`,1)
  post.rep.bc<-sapply(x,`[`,2)
  post.rep.b<-unname(sapply(post.rep.bc,`[`,1))
  post.rep.c<-unname(sapply(post.rep.bc,`[`,2))
  post.scr.m1.1.b<-sapply(post.screen.bc,`[`,1)
  post.scr.m1.b<-sapply(post.screen.bc,`[`,2)
  post.scr.m2.1.b<-sapply(post.screen.bc,`[`,3)
  post.scr.m2.b<-sapply(post.screen.bc,`[`,4)
  post.scr.msm1.b<-sapply(post.screen.bc,`[`,5)
  post.scr.msm2.b<-sapply(post.screen.bc,`[`,6)
  post.scr.f1.1.b<-sapply(post.screen.bc,`[`,7)
  post.scr.f1.b<-sapply(post.screen.bc,`[`,8)
  post.scr.f2.1.b<-sapply(post.screen.bc,`[`,9)
  post.scr.f2.b<-sapply(post.screen.bc,`[`,10)
  post.scr.m1.1.c<-sapply(post.screen.bc,`[`,11)
  post.scr.m1.c<-sapply(post.screen.bc,`[`,12)
  post.scr.m2.1.c<-sapply(post.screen.bc,`[`,13)
  post.scr.m2.c<-sapply(post.screen.bc,`[`,14)
  post.scr.msm1.c<-sapply(post.screen.bc,`[`,15)
  post.scr.msm2.c<-sapply(post.screen.bc,`[`,16)
  post.scr.f1.1.c<-sapply(post.screen.bc,`[`,17)
  post.scr.f1.c<-sapply(post.screen.bc,`[`,18)
  post.scr.f2.1.c<-sapply(post.screen.bc,`[`,19)
  post.scr.f2.c<-sapply(post.screen.bc,`[`,20)

  #rep.symp - reporting probability -  asymptomatic cases
  bez.rep.symp <- as.data.frame(cbind(bezA=ilogit(post.sample$theta[,"logit.rep.symp.a"]), bezD=ilogit(post.sample$theta[,"logit.rep.symp.d"]),bezB=post.rep.b, bezC=post.rep.c))
  rep.symp.post <-apply(bez.rep.symp, 1, bezier_fun)
  if (only_increasing_screening())
    rep.symp.post <- apply(rep.symp.post, 2, make_monotonically_increasing)
  rep.symp.post <- reshape2::melt(rep.symp.post)
  rep.symp.post$Var1 <- rep.symp.post$Var1+start.year-1
  rep.symp.prior.bez<- as.data.frame(cbind(a=rbeta(1000,prior.param1["rep.symp.a"],prior.param2["rep.symp.a"]),d=rbeta(1000,prior.param1["rep.symp.d"],prior.param2["rep.symp.d"]),b=rbeta(1000,prior.param1["rand.rep.symp.b"],prior.param2["rand.rep.symp.b"]),c=rbeta(1000,prior.param1["rand.rep.symp.c"],prior.param2["rand.rep.symp.c"])))
  rep.symp.prior.bez<-matrix(apply(rep.symp.prior.bez,1, prior_ctrl),ncol=4, byrow=TRUE)
  rep.symp.prior <-melt(apply(rep.symp.prior.bez, 1, bezier_fun))
  x<-data.frame(c(aggregate(value~Var1, rep.symp.prior, mean),aggregate(value~Var1,rep.symp.prior, min),  aggregate(value~Var1, rep.symp.prior, max)))
  x<-x[,c(1,2,4,6)]
  names(x)<-c("time", "mean", "min", "max")

  plot.rep.symp <- ggplot(data=rep.symp.post) +
    geom_line(aes(x=Var1, y=value, group=Var2), color="thistle4", alpha=0.5)+
    geom_ribbon(data=x, aes(x=seq(start.year, end.year),ymin=min, ymax=max), alpha=0.2)+
    stat_summary(fun.y = mean, geom='line', color = 'black', aes(x = Var1, y=value, group=1)) +
    theme_classic()+
    ylim(c(0,1)) +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8)) +
    labs(x="Year", y="Reporting probability")

  #transmission RR in MSM - change in condom use/behaviour in MSM leading to changes in transmission over time
  bez.behav <- as.data.frame(ilogit(post.sample$theta[,"logit.behav.lin"]))
  bez.behav <- apply(bez.behav, 1, behav_fun)
  for (i in 1:ncol(bez.behav)) {
    bez.behav[,i] <-
      ilogit(post.sample$theta[i,'logit.b.msm']) + (1-ilogit(post.sample$theta[i,'logit.b.msm']))*bez.behav[,i]
  }
  behav.post <-melt(bez.behav)
  behav.post$Var1 <- behav.post$Var1+start.year-1
  if (exists('gc_testing', envir=.GlobalEnv) && 'DelayedMSMTransmissionIncrease' %in% gc_testing) {
    behav.post$Var1 <- behav.post$Var1 + 8
  }
  # behav.post$value <- ilogit(post.sample$theta[,'logit.b.msm']) + (1-ilogit(post.sample$theta[,'logit.b.msm']))*behav.post$value
  # behav.prior.bez<- as.data.frame(rbeta(1000,prior.param1["behav.lin"],prior.param2["behav.lin"]))
  # behav.prior <-melt(apply(behav.prior.bez, 1, behav_fun))
  # x<-data.frame(c(aggregate(value~Var1, behav.prior, mean),aggregate(value~Var1,behav.prior, min),  aggregate(value~Var1, behav.prior, max)))
  # x<-x[,c(1,2,4,6)]
  # names(x)<-c("time", "mean", "min", "max")

  # plot.behav <- ggplot(data=behav.post) +
  #   geom_line(aes(x=Var1, y=value, group=Var2), color="lightcoral", alpha=0.5)+
  #   stat_summary(fun.y = mean, geom='line', color = 'black', aes(x = Var1, y=value, group=1), alpha=0.5) +
  #   # geom_ribbon(data=x, aes(x=seq(start.year, end.year), ymin=min, ymax=max), alpha=0.2)+
  #   theme_classic()+
  #   theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8)) +
  #   scale_y_continuous(breaks = seq(0, 1, 0.5)) +
  #   expand_limits(y=c(0,1)) +
  #   labs(x="Transmission Probability in MSM", y="probability")

  orig.behav.prior <- as.data.frame(cbind(
    x = seq(0, 1, 0.001),
    y = dbeta(seq(0, 1, 0.001), prior.param1["behav.lin"], prior.param2["behav.lin"])
  ))
  plot.orig.behav <-
    ggplot(
      data = data.frame(x=ilogit(post.sample$theta[,"logit.behav.lin"]))) +
    geom_area(
      data = orig.behav.prior,
      aes(x = x, y = y/max(y)),
      color = 'grey',
      alpha = 0.3
    ) +
    geom_histogram(aes(x=x, y=..ncount..), bins = 30, fill = 'lightcoral', color="lightcoral", alpha=0.6) +
    theme_classic() +
    theme(legend.position="none",
          axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          title=element_text(size=8),
          axis.ticks.x = element_blank()) +
    ylab("density") +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    xlab("Increase in Transmission in MSM\nRelative to (1-Transmission Pr)")

  # plot.rel.behav.msm <-
  #   ggplot(
  #     data = data.frame(x = ilogit(post.sample$theta[,'logit.behav.lin'])*(1-ilogit(post.sample$theta[,'logit.b.msm']))/ilogit(post.sample$theta[,'logit.b.msm']))) +
  #   geom_histogram(aes(x=x, y=..ncount..), bins = 30, fill = 'lightcoral', color="lightcoral", alpha=0.6) +
  #   theme_classic() +
  #   theme(legend.position="none",
  #         axis.text.x=element_text(size=8),
  #         axis.text.y=element_text(size=8),
  #         title=element_text(size=8),
  #         axis.ticks.x = element_blank()) +
  #   ylab("density") +
  #   xlab("Relative Increase in Transmission in MSM")

  # plot.abs.behav.msm <-
  #   ggplot(
  #     data = data.frame(x = ilogit(post.sample$theta[,'logit.behav.hetero'])*(1-ilogit(post.sample$theta[,'logit.b.msm'])))) +
  #   geom_histogram(aes(x=x, y=..ncount..), bins = 30, fill = 'lightcoral', color="lightcoral", alpha=0.6) +
  #   theme_classic() +
  #   theme(legend.position="none",
  #         axis.text.x=element_text(size=8),
  #         axis.text.y=element_text(size=8),
  #         title=element_text(size=8),
  #         axis.ticks.x = element_blank()) +
  #   scale_y_continuous(breaks = seq(0, 1, 0.5)) +
  #   ylab("density") +
  #   xlim(c(0,1)) +
    # xlab("Transmission Increase in M to M")


  b.msm.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["b.msm"],prior.param2["b.msm"])))

  b.inc.msm.prior <- data.frame(
    x = rbeta(10000, shape1 = prior.param1['b.msm'], shape2 = prior.param2['b.msm']),
    y = rbeta(10000, shape1 = prior.param1['behav.lin'], shape2 = prior.param2['behav.lin'])
  )

  b.inc.msm.prior %<>% mutate(
    z = x + (1-x)*y
  )

  b.inc.msm.prior <- MASS::fitdistr(b.inc.msm.prior$z, densfun = 'beta', start = list(shape1=prior.param1['b.msm'], shape2=prior.param2['b.msm']))

  b.inc.msm.prior <- data.frame(x=seq(0,1,0.01), y = dbeta(seq(0,1,0.01), b.inc.msm.prior$estimate[[1]], b.inc.msm.prior$estimate[[2]]))

  plot.behav.msm.inc <-
    ggplot(
      data = data.frame(
        x =  ilogit(post.sample$theta[,'logit.b.msm']) + ilogit(post.sample$theta[,'logit.behav.lin'])*(1-ilogit(post.sample$theta[,'logit.b.msm'])))) +
    geom_area(data = b.inc.msm.prior, aes(x=x, y = y/max(y)), color = 'grey', fill = 'grey', alpha=0.6, size=0) +
    geom_histogram(aes(x=x, y=..ncount..), fill="dodgerblue3",size=0.5, colour="dodgerblue3", alpha=0.6, binwidth=0.025)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,1))+
    labs(x=paste0("Transmission probability in ", gc_env$end.year, ":\n M to M"))


  b.inc.m.prior <- data.frame(
    x = rbeta(2000, shape1 = prior.param1['b.m'], shape2 = prior.param2['b.m']),
    y = rbeta(2000, shape1 = prior.param1['behav.hetero'], shape2 = prior.param2['behav.hetero'])
  )

  b.inc.m.prior %<>% mutate(
    z = x + (1-x)*y
  ) %>% filter(is.finite(z))

  b.inc.m.prior <- MASS::fitdistr(b.inc.m.prior$z, densfun = 'beta', start = list(shape1=prior.param1['b.m'], shape2=prior.param2['b.m']))

  b.inc.m.prior <- data.frame(x=seq(0,1,0.01), y = dbeta(seq(0,1,0.01), b.inc.m.prior$estimate[[1]], b.inc.m.prior$estimate[[2]]))


  plot.behav.m.inc <-
    ggplot(
      data = data.frame(
        x =  ilogit(post.sample$theta[,'logit.b.m']) + ilogit(post.sample$theta[,'logit.behav.hetero'])*(1-ilogit(post.sample$theta[,'logit.b.m'])))) +
    geom_area(data = b.inc.m.prior, aes(x=x, y = y/max(y)), color = 'grey', fill = 'grey', alpha=0.6, size=0) +
    geom_histogram(aes(x=x, y=..ncount..), fill="dodgerblue3",size=0.5, colour="dodgerblue3", alpha=0.6, binwidth=0.025)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,1))+
    labs(x=paste0("Transmission probability ", gc_env$end.year, ":\n F to M"))

  b.inc.f.prior <- data.frame(
    x = rbeta(2000, shape1 = prior.param1['b.f'], shape2 = prior.param2['b.f']),
    y = rbeta(2000, shape1 = prior.param1['behav.hetero'], shape2 = prior.param2['behav.hetero'])
  )

  b.inc.f.prior %<>% mutate(
    z = x + (1-x)*y
  )

  b.inc.f.prior <- MASS::fitdistr(b.inc.f.prior$z, densfun = 'beta', start = list(shape1=prior.param1['b.f'], shape2=prior.param2['b.f']))

  b.inc.f.prior <- data.frame(x=seq(0,1,0.01), y = dbeta(seq(0,1,0.01), b.inc.f.prior$estimate[[1]], b.inc.f.prior$estimate[[2]]))


  plot.behav.f.inc <-
    ggplot(
      data = data.frame(
        x =  ilogit(post.sample$theta[,'logit.b.f']) + ilogit(post.sample$theta[,'logit.behav.hetero'])*(1-ilogit(post.sample$theta[,'logit.b.f'])))) +
    geom_area(data = b.inc.f.prior, aes(x=x, y = y/max(y)), color = 'grey', fill = 'grey', alpha=0.6, size=0) +
    geom_histogram(aes(x=x, y=..ncount..), fill="dodgerblue3",size=0.5, colour="dodgerblue3", alpha=0.6, binwidth=0.025)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=6), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    ylab("density") +
    expand_limits(x=c(0,1))+
    labs(x=paste0("Transmission probability ", gc_env$end.year, ":\n F to M"))

    #   ggplot(
    #     data = data.frame(x = ilogit(post.sample$theta[,'logit.behav.hetero'])*(1-ilogit(post.sample$theta[,'logit.b.msm'])))) +

  # plot.rel.behav.m <-
  #   ggplot(
  #     data = data.frame(x = ilogit(post.sample$theta[,'logit.behav.hetero'])*(1-ilogit(post.sample$theta[,'logit.b.m']))/ilogit(post.sample$theta[,'logit.b.m']))) +
  #   geom_histogram(aes(x=x, y=..ncount..), bins = 30, fill = 'lightcoral', color="lightcoral", alpha=0.6) +
  #   theme_classic() +
  #   theme(legend.position="none",
  #         axis.text.x=element_text(size=8),
  #         axis.text.y=element_text(size=8),
  #         title=element_text(size=8),
  #         axis.ticks.x = element_blank()) +
  #   ylab("density") +
  #   xlab("Relative Increase in Transmission in M to F")



  # plot.rel.behav.f <-
  #   ggplot(
  #     data = data.frame(x = ilogit(post.sample$theta[,'logit.behav.hetero'])*(1-ilogit(post.sample$theta[,'logit.b.f']))/ilogit(post.sample$theta[,'logit.b.f']))) +
  #   geom_histogram(aes(x=x, y=..ncount..), bins = 30, fill = 'lightcoral', color="lightcoral", alpha=0.6) +
  #   theme_classic() +
  #   theme(legend.position="none",
  #         axis.text.x=element_text(size=8),
  #         axis.text.y=element_text(size=8),
  #         title=element_text(size=8),
  #         axis.ticks.x = element_blank()) +
  #   ylab("density") +
  #   xlab("Relative Increase in Transmission in F to M")



  if (check_variation("separate_increased_risk")) {
    plot.behav.hetero <-
      ggplot(
        data = data.frame(x = ilogit(post.sample$theta[,'logit.behav.hetero']))) +
      geom_area(
        data = cbind.data.frame(
          x = seq(0,1,0.001),
          y = dbeta(seq(0,1,0.001), prior.param1['behav.hetero'], prior.param2['behav.hetero'])),
        aes(x = x, y = y/max(y)),
        color = 'grey',
        alpha = 0.3
      ) +
      geom_histogram(aes(x=x, y=..ncount..), bins = 30, fill = 'lightcoral', color="lightcoral", alpha=0.6) +
      theme_classic() +
      theme(legend.position="none",
            axis.text.x=element_text(size=8),
            axis.text.y=element_text(size=8),
            title=element_text(size=8),
            axis.ticks.x = element_blank()) +
      ylab("density") +
      scale_y_continuous(breaks = seq(0, 1, 0.5)) +
      xlab("Increase in Transmission in MSW and W\nRelative to (1-Transmission Pr)")

    # plot.abs.behav.m <-
    #   ggplot(
    #     data = data.frame(x = ilogit(post.sample$theta[,'logit.behav.hetero'])*(1-ilogit(post.sample$theta[,'logit.b.m'])))) +
    #   geom_histogram(aes(x=x, y=..ncount..), bins = 30, fill = 'lightcoral', color="lightcoral", alpha=0.6) +
    #   theme_classic() +
    #   theme(legend.position="none",
    #         axis.text.x=element_text(size=8),
    #         axis.text.y=element_text(size=8),
    #         title=element_text(size=8),
    #         axis.ticks.x = element_blank()) +
    #   scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    #   ylab("density") +
    #   xlim(c(0,1)) +
    #   xlab("Transmission Increase in F to M")

    # plot.abs.behav.f <-
    #   ggplot(
    #     data = data.frame(x = ilogit(post.sample$theta[,'logit.behav.hetero'])*(1-ilogit(post.sample$theta[,'logit.b.f'])))) +
    #   geom_histogram(aes(x=x, y=..ncount..), bins = 30, fill = 'lightcoral', color="lightcoral", alpha=0.6) +
    #   theme_classic() +
    #   theme(legend.position="none",
    #         axis.text.x=element_text(size=8),
    #         axis.text.y=element_text(size=8),
    #         title=element_text(size=8),
    #         axis.ticks.x = element_blank()) +
    #   scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    #   ylab("density") +
    #   xlim(c(0,1)) +
    #   xlab("Transmission Increase in M to F")


  } else {
    plot.behav.hetero <- blankPlot
    plot.abs.behav.m <- blankPlot
    plot.abs.behav.f <- blankPlot
  }


  # we're going to use this named_numeric_from_list function inside
  # the construction of the 3d screening array.
  named_numeric_from_list <- function(l) {
    n <- names(l)
    l <- as.numeric(l)
    names(l) <- n
    return(l)
  }

  # construct the trace as a data frame
  trace <- as.data.frame(post.sample$theta.list)

  # this is a list of screening arrays.
  screening_arrays <-
    sapply(1:nrow(trace), function(i)
      construct_screening_matrix(named_numeric_from_list(trace[i, ]), update_ctrl(named_numeric_from_list(trace[i, ]))),
      simplify = FALSE)

  # turn the list of arrays into a 3d array
  # it's dimensions are:
  # 1 - time
  # 2 - population
  # 3 - simulation
  screening <- abind::abind(screening_arrays, along = 3)

  # switch
  screening <- aperm(screening, c(2, 1, 3))

  #screen.f1 - screening rate in youngest females (other/Hispanic)
  # bez.screen.f1 <- as.data.frame(cbind(bezA=ilogit(post.sample$theta[,"logit.screen.f1.a"]), bezD=ilogit(post.sample$theta[,"logit.screen.f1.d"]),bezB= post.scr.f1.b, bezC=post.scr.f1.c))
  # screen.f1.post <-melt(apply(bez.screen.f1, 1, bezier_fun))
  screen.f1.post <- screening['female.other.young.low', ,]
  screen.f1.post <- melt(screen.f1.post)
  screen.f1.post$Var1 <- screen.f1.post$Var1+start.year-1
  screen.f1.prior.bez<- as.data.frame(cbind(a=rbeta(1000,prior.param1["screen.f1.a"],prior.param2["screen.f1.a"]),d=rbeta(1000,prior.param1["screen.f1.d"],prior.param2["screen.f1.d"]),b=rbeta(1000,prior.param1["rand.screen.f1.b"],prior.param2["rand.screen.f1.b"]),c=rbeta(1000,prior.param1["rand.screen.f1.c"],prior.param2["rand.screen.f1.c"])))
  screen.f1.prior.bez <- matrix(apply(screen.f1.prior.bez,1, prior_ctrl),ncol=4, byrow=TRUE)
  screen.f1.prior <-melt(apply(screen.f1.prior.bez, 1, bezier_fun))
  x<-data.frame(c(aggregate(value~Var1, screen.f1.prior, mean),aggregate(value~Var1,screen.f1.prior, min),  aggregate(value~Var1, screen.f1.prior, max)))
  x<-x[,c(1,2,4,6)]
  names(x)<-c("time", "mean", "min", "max")

  plot.screen.f1 <- ggplot(data=screen.f1.post) +
    geom_line(aes(x=Var1, y=value, group=Var2), color="steelblue1", alpha=0.5)+
    geom_ribbon(data=x, aes(x=seq(start.year, end.year,1),ymin=min, ymax=max), alpha=0.2)+
    stat_summary(fun.y = mean, geom='line', color = 'black', aes(x = Var1, y=value, group=1), alpha=0.5) +
    theme_classic()+
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8)) +
    coord_cartesian(ylim=c(0,0.8))+
    labs(title="Other and Hispanic F:\n 15-24y", x="Year", y="Asymptomatic screen\nand treat rate")

  #screen.f1.1 - screening rate in youngest females (black)
  # bez.screen.f1.1 <- as.data.frame(cbind(bezA=ilogit(post.sample$theta[,"logit.screen.f1.1.a"]), bezD=ilogit(post.sample$theta[,"logit.screen.f1.1.d"]),bezB= post.scr.f1.1.b, bezC=post.scr.f1.1.c))
  # screen.f1.1.post <-melt(apply(bez.screen.f1.1, 1, bezier_fun))
  screen.f1.1.post <- screening['female.black.young.low', ,]
  screen.f1.1.post <- melt(screen.f1.1.post)
  screen.f1.1.post$Var1 <- screen.f1.1.post$Var1+start.year-1
  screen.f1.1.prior.bez<- as.data.frame(cbind(a=rbeta(1000,prior.param1["screen.f1.1.a"],prior.param2["screen.f1.1.a"]),d=rbeta(1000,prior.param1["screen.f1.1.d"],prior.param2["screen.f1.1.d"]),b=rbeta(1000,prior.param1["rand.screen.f1.1.b"],prior.param2["rand.screen.f1.1.b"]),c=rbeta(1000,prior.param1["rand.screen.f1.1.c"],prior.param2["rand.screen.f1.1.c"])))
  screen.f1.1.prior.bez <- matrix(apply(screen.f1.1.prior.bez,1, prior_ctrl),ncol=4, byrow=TRUE)
  screen.f1.1.prior <-melt(apply(screen.f1.1.prior.bez, 1, bezier_fun))
  x<-data.frame(c(aggregate(value~Var1, screen.f1.1.prior, mean),aggregate(value~Var1,screen.f1.1.prior, min),  aggregate(value~Var1, screen.f1.1.prior, max)))
  x<-x[,c(1,2,4,6)]
  names(x)<-c("time", "mean", "min", "max")

  plot.screen.f1.1 <- ggplot(data=screen.f1.1.post) +
    geom_line(aes(x=Var1, y=value, group=Var2), color="steelblue1", alpha=0.5)+
    geom_ribbon(data=x, aes(x=seq(start.year, end.year,1),ymin=min, ymax=max), alpha=0.2)+
    stat_summary(fun.y = mean, geom='line', color = 'black', aes(x = Var1, y=value, group=1), alpha=0.5) +
    theme_classic()+
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8)) +
    coord_cartesian(ylim=c(0,0.8))+
    labs(title="Black F:\n 15-24y", x="Year", y="Asymptomatic screen\nand treat rate")

  #screen.f2 - screening rate in older females (Hispanic and other)
  # bez.screen.f2 <- as.data.frame(cbind(bezA=ilogit(post.sample$theta[,"logit.screen.f2.a"]), bezD=ilogit(post.sample$theta[,"logit.screen.f2.d"]),bezB=post.scr.f2.b, bezC=post.scr.f2.c))
  # screen.f2.post <-melt(apply(bez.screen.f2, 1, bezier_fun))
  screen.f2.post <- screening['female.other.old.low', , ]
  screen.f2.post <- melt(screen.f2.post)
  screen.f2.post$Var1 <- screen.f2.post$Var1+start.year-1
  screen.f2.prior.bez<- as.data.frame(cbind(a=rbeta(1000,prior.param1["screen.f2.a"],prior.param2["screen.f2.a"]),d=rbeta(1000,prior.param1["screen.f2.d"],prior.param2["screen.f2.d"]),b=rbeta(1000,prior.param1["rand.screen.f2.b"],prior.param2["rand.screen.f2.b"]),c=rbeta(1000,prior.param1["rand.screen.f2.c"],prior.param2["rand.screen.f2.c"])))
  screen.f2.prior.bez <- matrix(apply(screen.f2.prior.bez,1, prior_ctrl),ncol=4, byrow=TRUE)
  screen.f2.prior <-melt(apply(screen.f2.prior.bez, 1, bezier_fun))
  x<-data.frame(c(aggregate(value~Var1, screen.f2.prior, mean),aggregate(value~Var1,screen.f2.prior, min),  aggregate(value~Var1, screen.f2.prior, max)))
  x<-x[,c(1,2,4,6)]
  names(x)<-c("time", "mean", "min", "max")

  plot.screen.f2 <- ggplot(data=screen.f2.post) +
    geom_line(aes(x=Var1, y=value, group=Var2), color="steelblue1", alpha=0.5)+
    geom_ribbon(data=x, aes(x=seq(start.year, end.year,1),ymin=min, ymax=max), alpha=0.2)+
    stat_summary(fun.y = mean, geom='line', color = 'black', aes(x = Var1, y=value, group=1)) +
    theme_classic()+
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8)) +
    coord_cartesian(ylim=c(0,0.8))+
    labs(title="Other and Hispanic F:\n 25-39y",x="Year", y="Asymptomatic screen\nand treat rate")

  #screen.f2.1 - screening rate in older females (black)
  # bez.screen.f2.1 <- as.data.frame(cbind(bezA=ilogit(post.sample$theta[,"logit.screen.f2.1.a"]), bezD=ilogit(post.sample$theta[,"logit.screen.f2.1.d"]),bezB=post.scr.f2.1.b, bezC=post.scr.f2.1.c))
  # screen.f2.1.post <-melt(apply(bez.screen.f2.1, 1, bezier_fun))
  screen.f2.1.post <- screening['female.black.old.low', , ]
  screen.f2.1.post <- melt(screen.f2.1.post)
  screen.f2.1.post$Var1 <- screen.f2.1.post$Var1+start.year-1
  screen.f2.1.prior.bez<- as.data.frame(cbind(a=rbeta(1000,prior.param1["screen.f2.1.a"],prior.param2["screen.f2.1.a"]),d=rbeta(1000,prior.param1["screen.f2.1.d"],prior.param2["screen.f2.1.d"]),b=rbeta(1000,prior.param1["rand.screen.f2.1.b"],prior.param2["rand.screen.f2.1.b"]),c=rbeta(1000,prior.param1["rand.screen.f2.1.c"],prior.param2["rand.screen.f2.1.c"])))
  screen.f2.1.prior.bez <- matrix(apply(screen.f2.1.prior.bez,1, prior_ctrl),ncol=4, byrow=TRUE)
  screen.f2.1.prior <-melt(apply(screen.f2.1.prior.bez, 1, bezier_fun))
  x<-data.frame(c(aggregate(value~Var1, screen.f2.1.prior, mean),aggregate(value~Var1,screen.f2.1.prior, min),  aggregate(value~Var1, screen.f2.1.prior, max)))
  x<-x[,c(1,2,4,6)]
  names(x)<-c("time", "mean", "min", "max")

  plot.screen.f2.1 <- ggplot(data=screen.f2.1.post) +
    geom_line(aes(x=Var1, y=value, group=Var2), color="steelblue1", alpha=0.5)+
    geom_ribbon(data=x, aes(x=seq(start.year, end.year,1),ymin=min, ymax=max), alpha=0.2)+
    stat_summary(fun.y = mean, geom='line', color = 'black', aes(x = Var1, y=value, group=1), alpha=0.5) +
    theme_classic()+
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8)) +
    coord_cartesian(ylim=c(0,0.8))+
    labs(title="Black F:\n 25-39y",x="Year", y="Asymptomatic screen\nand treat rate")

  #screen.m1 - screening rate in youngest males (other/Hispanic)
  # bez.screen.m1 <- as.data.frame(cbind(bezA=ilogit(post.sample$theta[,"logit.screen.m1.a"]), bezD=ilogit(post.sample$theta[,"logit.screen.m1.d"]),bezB= post.scr.m1.b, bezC=post.scr.m1.c))
  # screen.m1.post <-melt(apply(bez.screen.m1, 1, bezier_fun))
  screen.m1.post <- screening['male.other.young.low', ,]
  screen.m1.post <- melt(screen.m1.post)
  screen.m1.post$Var1 <- screen.m1.post$Var1+start.year-1
  screen.m1.prior.bez<- as.data.frame(cbind(a=rbeta(1000,prior.param1["screen.m1.a"],prior.param2["screen.m1.a"]),d=rbeta(1000,prior.param1["screen.m1.d"],prior.param2["screen.m1.d"]),b=rbeta(1000,prior.param1["rand.screen.m1.b"],prior.param2["rand.screen.m1.b"]),c=rbeta(1000,prior.param1["rand.screen.m1.c"],prior.param2["rand.screen.m1.c"])))
  screen.m1.prior.bez <- matrix(apply(screen.m1.prior.bez,1, prior_ctrl),ncol=4, byrow=TRUE)
  screen.m1.prior <-melt(apply(screen.m1.prior.bez, 1, bezier_fun))
  x<-data.frame(c(aggregate(value~Var1, screen.m1.prior, mean),aggregate(value~Var1,screen.m1.prior, min),  aggregate(value~Var1, screen.m1.prior, max)))
  x<-x[,c(1,2,4,6)]
  names(x)<-c("time", "mean", "min", "max")

  plot.screen.m1 <- ggplot(data=screen.m1.post) +
    geom_line(aes(x=Var1, y=value, group=Var2), color="steelblue1", alpha=0.5)+
    geom_ribbon(data=x, aes(x=seq(start.year, end.year,1),ymin=min, ymax=max), alpha=0.2)+
    stat_summary(fun.y = mean, geom='line', color = 'black', aes(x = Var1, y=value, group=1), alpha=0.5) +
    theme_classic()+
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8)) +
    coord_cartesian(ylim=c(0,0.8))+
    labs(title="Other and Hispanic M:\n 15-24y",x="Year", y="Asymptomatic screen\nand treat rate")

  #screen.m2 - screening rate in older males (other/Hispanic)
  # bez.screen.m2 <- as.data.frame(cbind(bezA=ilogit(post.sample$theta[,"logit.screen.m2.a"]), bezD=ilogit(post.sample$theta[,"logit.screen.m2.d"]),bezB=post.scr.m2.b, bezC=post.scr.m2.c))
  # screen.m2.post <-melt(apply(bez.screen.m2, 1, bezier_fun))
  screen.m2.post <- screening['male.other.old.low', ,]
  screen.m2.post <- melt(screen.m2.post)
  screen.m2.post$Var1 <- screen.m2.post$Var1+start.year-1
  screen.m2.prior.bez<- as.data.frame(cbind(a=rbeta(1000,prior.param1["screen.m2.a"],prior.param2["screen.m2.a"]),d=rbeta(1000,prior.param1["screen.m2.d"],prior.param2["screen.m2.d"]),b=rbeta(1000,prior.param1["rand.screen.m2.b"],prior.param2["rand.screen.m2.b"]),c=rbeta(1000,prior.param1["rand.screen.m2.c"],prior.param2["rand.screen.m2.c"])))
  screen.m2.prior.bez <- matrix(apply(screen.m2.prior.bez,1, prior_ctrl),ncol=4, byrow=TRUE)
  screen.m2.prior <-melt(apply(screen.m2.prior.bez, 1, bezier_fun))
  x<-data.frame(c(aggregate(value~Var1, screen.m2.prior, mean),aggregate(value~Var1,screen.m2.prior, min),  aggregate(value~Var1, screen.m2.prior, max)))
  x<-x[,c(1,2,4,6)]
  names(x)<-c("time", "mean", "min", "max")

  plot.screen.m2 <- ggplot(data=screen.m2.post) +
    geom_line(aes(x=Var1, y=value, group=Var2), color="steelblue1", alpha=0.5)+
    geom_ribbon(data=x, aes(x=seq(start.year, end.year,1),ymin=min, ymax=max), alpha=0.2)+
    stat_summary(fun.y = mean, geom='line', color = 'black', aes(x = Var1, y=value, group=1), alpha=0.5) +
    theme_classic()+
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8)) +
    coord_cartesian(ylim=c(0,0.8))+
    labs(title="Other and Hispanic M:\n 25-39y", x="Year", y="Asymptomatic screen\nand treat rate")

  #screen.m1.1 - screening rate in youngest males (black)
  # bez.screen.m1.1 <- as.data.frame(cbind(bezA=ilogit(post.sample$theta[,"logit.screen.m1.1.a"]), bezD=ilogit(post.sample$theta[,"logit.screen.m1.1.d"]),bezB= post.scr.m1.1.b, bezC=post.scr.m1.1.c))
  # screen.m1.1.post <-melt(apply(bez.screen.m1.1, 1, bezier_fun))
  screen.m1.1.post <- screening['male.black.young.low', , ]
  screen.m1.1.post <- melt(screen.m1.1.post)
  screen.m1.1.post$Var1 <- screen.m1.1.post$Var1+start.year-1
  screen.m1.1.prior.bez<- as.data.frame(cbind(a=rbeta(1000,prior.param1["screen.m1.1.a"],prior.param2["screen.m1.1.a"]),d=rbeta(1000,prior.param1["screen.m1.1.d"],prior.param2["screen.m1.1.d"]),b=rbeta(1000,prior.param1["rand.screen.m1.1.b"],prior.param2["rand.screen.m1.1.b"]),c=rbeta(1000,prior.param1["rand.screen.m1.1.c"],prior.param2["rand.screen.m1.1.c"])))
  screen.m1.1.prior.bez <- matrix(apply(screen.m1.1.prior.bez,1, prior_ctrl),ncol=4, byrow=TRUE)
  screen.m1.1.prior <-melt(apply(screen.m1.1.prior.bez, 1, bezier_fun))
  x<-data.frame(c(aggregate(value~Var1, screen.m1.1.prior, mean),aggregate(value~Var1,screen.m1.1.prior, min),  aggregate(value~Var1, screen.m1.1.prior, max)))
  x<-x[,c(1,2,4,6)]
  names(x)<-c("time", "mean", "min", "max")

  plot.screen.m1.1 <- ggplot(data=screen.m1.1.post) +
    geom_line(aes(x=Var1, y=value, group=Var2), color="steelblue1",alpha=0.5)+
    geom_ribbon(data=x, aes(x=time+start.year-1,ymin=min, ymax=max), alpha=0.2)+
    stat_summary(fun.y = mean, geom='line', color = 'black', aes(x = Var1, y=value, group=1), alpha=0.5) +
    theme_classic()+
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8)) +
    coord_cartesian(ylim=c(0,0.8))+
    labs(title="Black M:\n 15-24y",x="Year", y="Asymptomatic screen\nand treat rate")

  #screen.m2.1 - screening rate in older males (black)
  # bez.screen.m2.1 <- as.data.frame(cbind(bezA=ilogit(post.sample$theta[,"logit.screen.m2.1.a"]), bezD=ilogit(post.sample$theta[,"logit.screen.m2.1.d"]),bezB=post.scr.m2.1.b, bezC=post.scr.m2.1.c))
  # screen.m2.1.post <-melt(apply(bez.screen.m2.1, 1, bezier_fun))
  screen.m2.1.post <- screening['male.black.old.low', ,]
  screen.m2.1.post <- melt(screen.m2.1.post)
  screen.m2.1.post$Var1 <- screen.m2.1.post$Var1+start.year-1
  screen.m2.1.prior.bez<- as.data.frame(cbind(a=rbeta(1000,prior.param1["screen.m2.1.a"],prior.param2["screen.m2.1.a"]),d=rbeta(1000,prior.param1["screen.m2.1.d"],prior.param2["screen.m2.1.d"]),b=rbeta(1000,prior.param1["rand.screen.m2.1.b"],prior.param2["rand.screen.m2.1.b"]),c=rbeta(1000,prior.param1["rand.screen.m2.1.c"],prior.param2["rand.screen.m2.1.c"])))
  screen.m2.1.prior.bez <- matrix(apply(screen.m2.1.prior.bez,1, prior_ctrl),ncol=4, byrow=TRUE)
  screen.m2.1.prior <-melt(apply(screen.m2.1.prior.bez, 1, bezier_fun))
  x<-data.frame(c(aggregate(value~Var1, screen.m2.1.prior, mean),aggregate(value~Var1,screen.m2.1.prior, min),  aggregate(value~Var1, screen.m2.1.prior, max)))
  x<-x[,c(1,2,4,6)]
  names(x)<-c("time", "mean", "min", "max")

  plot.screen.m2.1 <- ggplot(data=screen.m2.1.post) +
    geom_line(aes(x=Var1, y=value, group=Var2), color="steelblue1",alpha=0.5)+
    geom_ribbon(data=x, aes(x=seq(start.year, end.year,1),ymin=min, ymax=max), alpha=0.2)+
    stat_summary(fun.y = mean, geom='line', color = 'black', aes(x = Var1, y=value, group=1), alpha=0.5) +
    theme_classic()+
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8)) +
    coord_cartesian(ylim=c(0,0.8))+
    labs(title="Black M:\n 25-39y", x="Year", y="Asymptomatic screen\nand treat rate")

  #screen.msm1 - screening rate in youngest MSM
  # bez.screen.msm1 <- as.data.frame(cbind(bezA=ilogit(post.sample$theta[,"logit.screen.msm1.a"]), bezD=ilogit(post.sample$theta[,"logit.screen.msm1.d"]),bezB= post.scr.msm1.b, bezC=post.scr.msm1.c))
  # if (exists('gc_testing', envir = .GlobalEnv) && 'flat_msm_screening' %in% gc_testing) {
  # bez.screen.msm1 <- as.data.frame(cbind(bezA=ilogit(post.sample$theta[,"logit.screen.msm1.a"]), bezD=ilogit(post.sample$theta[,"logit.screen.msm1.a"]),bezB= ilogit(post.sample$theta[,"logit.screen.msm1.a"]), bezC=ilogit(post.sample$theta[,"logit.screen.msm1.a"])))
  # }
  # screen.msm1.post <-melt(apply(bez.screen.msm1, 1, bezier_fun))
  screen.msm1.post <- screening['male.msm.young.low', , ]
  screen.msm1.post <- melt(screen.msm1.post)
  screen.msm1.post$Var1 <- screen.msm1.post$Var1+start.year-1
  screen.msm1.prior.bez<- as.data.frame(cbind(a=rbeta(1000,prior.param1["screen.msm1.a"],prior.param2["screen.msm1.a"]),d=rbeta(1000,prior.param1["screen.msm1.d"],prior.param2["screen.msm1.d"]),b=rbeta(1000,prior.param1["rand.screen.msm1.b"],prior.param2["rand.screen.msm1.b"]),c=rbeta(1000,prior.param1["rand.screen.msm1.c"],prior.param2["rand.screen.msm1.c"])))
  screen.msm1.prior.bez <- matrix(apply(screen.msm1.prior.bez,1, prior_ctrl),ncol=4, byrow=TRUE)
  screen.msm1.prior <-melt(apply(screen.msm1.prior.bez, 1, bezier_fun))
  x<-data.frame(c(aggregate(value~Var1, screen.msm1.prior, mean),aggregate(value~Var1,screen.msm1.prior, min),  aggregate(value~Var1, screen.msm1.prior, max)))
  x<-x[,c(1,2,4,6)]
  names(x)<-c("time", "mean", "min", "max")

  plot.screen.msm1 <- ggplot(data=screen.msm1.post) +
    geom_line(aes(x=Var1, y=value, group=Var2), color="steelblue1", alpha=0.5)+
    geom_ribbon(data=x, aes(x=seq(start.year, end.year,1),ymin=min, ymax=max), alpha=0.2)+
    stat_summary(fun.y = mean, geom='line', color = 'black', aes(x = Var1, y=value, group=1), alpha=0.5) +
    theme_classic()+
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8)) +
    coord_cartesian(ylim=c(0,0.8))+
    labs(title="MSM:\n 15-24y", x="Year", y="Asymptomatic screen\nand treat rate")

  #screen.msm2 - screening rate in older MSM
  # bez.screen.msm2 <- as.data.frame(cbind(bezA=ilogit(post.sample$theta[,"logit.screen.msm2.a"]), bezD=ilogit(post.sample$theta[,"logit.screen.msm2.d"]),bezB=post.scr.msm2.b, bezC=post.scr.msm2.c))
  # if (exists('gc_testing', envir = .GlobalEnv) && 'flat_msm_screening' %in% gc_testing) {
  #   bez.screen.msm2 <- as.data.frame(cbind(bezA=ilogit(post.sample$theta[,"logit.screen.msm2.a"]), bezD=ilogit(post.sample$theta[,"logit.screen.msm2.a"]),bezB= ilogit(post.sample$theta[,"logit.screen.msm2.a"]), bezC=ilogit(post.sample$theta[,"logit.screen.msm2.a"])))
  # }
  # screen.msm2.post <-melt(apply(bez.screen.msm2, 1, bezier_fun))
  screen.msm2.post <- screening['male.msm.old.low', , ]
  screen.msm2.post <- melt(screen.msm2.post)
  screen.msm2.post$Var1 <- screen.msm2.post$Var1+start.year-1
  screen.msm2.prior.bez<- as.data.frame(cbind(a=rbeta(1000,prior.param1["screen.msm2.a"],prior.param2["screen.msm2.a"]),d=rbeta(1000,prior.param1["screen.msm2.d"],prior.param2["screen.msm2.d"]),b=rbeta(1000,prior.param1["rand.screen.msm2.b"],prior.param2["rand.screen.msm2.b"]),c=rbeta(1000,prior.param1["rand.screen.msm2.c"],prior.param2["rand.screen.msm2.c"])))
  screen.msm2.prior.bez <- matrix(apply(screen.msm2.prior.bez,1, prior_ctrl),ncol=4, byrow=TRUE)
  screen.msm2.prior <-melt(apply(screen.msm2.prior.bez, 1, bezier_fun))
  x<-data.frame(c(aggregate(value~Var1, screen.msm2.prior, mean),aggregate(value~Var1,screen.msm2.prior, min),  aggregate(value~Var1, screen.msm2.prior, max)))
  x<-x[,c(1,2,4,6)]
  names(x)<-c("time", "mean", "min", "max")

  plot.screen.msm2 <- ggplot(data=screen.msm2.post) +
    geom_line(aes(x=Var1, y=value, group=Var2), color="steelblue1", alpha=0.5)+
    geom_ribbon(data=x, aes(x=seq(start.year, end.year,1),ymin=min, ymax=max), alpha=0.2)+
    stat_summary(fun.y = mean, geom='line', color = 'black', aes(x = Var1, y=value, group=1), alpha=0.5) +
    theme_classic()+
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8)) +
    coord_cartesian(ylim=c(0,0.8))+
    labs(title="MSM:\n 25-39y", x="Year", y="Asymptomatic screen\nand treat rate")

  ### save plots to pdf ###
  if (! dir.exists(output_directory))
    stop("directory specified does not exist.")

  pdf(file=file.path(output_directory, filename), width=11, height=8.5, paper="a4r")

  grid.arrange(plot.prev,     plot.age,     plot.subpop,
               plot.p.msm.y,  plot.p.msm.o, plot.p.symp, ncol=3)

  grid.arrange(plot.diag.y.m, plot.diag.y.m.1, plot.diag.y.m.2, plot.diag.y.m.3, plot.diag.y.msm,
               plot.diag.o.m, plot.diag.o.m.1, plot.diag.o.m.2, plot.diag.o.m.3, plot.diag.o.msm,
               plot.diag.y.f, plot.diag.y.f.1, plot.diag.y.f.2, plot.diag.y.f.3, blankPlot,
               plot.diag.o.f, plot.diag.o.f.1, plot.diag.o.f.2, plot.diag.o.f.3, blankPlot,  ncol=5)

  grid.arrange(plot.rr.diag.m.y.1, plot.rr.diag.f.y.1, plot.rr.diag.m.y.2, plot.rr.diag.f.y.2,
               plot.rr.diag.m.o.1, plot.rr.diag.f.o.1, plot.rr.diag.m.o.2, plot.rr.diag.f.o.2, ncol=4)

  grid.arrange(plot.inc.y.m.1, plot.inc.y.m.2, plot.inc.y.m.3, plot.inc.y.msm,
               plot.inc.o.m.1, plot.inc.o.m.2, plot.inc.o.m.3, plot.inc.o.msm,
               plot.inc.y.f.1, plot.inc.y.f.2, plot.inc.y.f.3, plot.inc.msm.mean,
               plot.inc.o.f.1, plot.inc.o.f.2, plot.inc.o.f.3, blankPlot, ncol=4)

  grid.arrange(plot.screen.m1.1,  plot.screen.m1,      plot.screen.msm1,   plot.screen.f1.1,    plot.screen.f1,
               plot.screen.m2.1,  plot.screen.m2,      plot.screen.msm2,   plot.screen.f2.1,    plot.screen.f2,
               plot.rr.screen.m3, plot.rr.screen.f3,   plot.rr.screen.ac,  blankPlot,           blankPlot,
               plot.rep.symp,     plot.rr.rep.symp.m1 ,plot.rr.rep.symp.m, plot.rr.rep.symp.f1, plot.rr.rep.symp.f,
               bottom=textGrob("Priors/posteriors (1/3)", x=0.01, y=1, just="left"),
               ncol=5)

  grid.arrange(plot.b.m,             plot.b.f,             plot.b.msm,             plot.behav.hetero,
               plot.behav.m.inc,     plot.behav.f.inc,     plot.behav.msm.inc,     plot.orig.behav,
               plot.symp.m,          plot.symp.f,          plot.symp.msm,          blankPlot,
               plot.dur.inf.symp.m,  plot.dur.inf.symp.f,  plot.dur.inf.symp.msm,  blankPlot,
               plot.dur.inf.asymp.m, plot.dur.inf.asymp.f, plot.dur.inf.asymp.msm, blankPlot,
               plot.theta.1,         plot.theta.2,         plot.theta.3,           plot.theta.4,
               plot.theta.5,         plot.theta.6,         plot.theta.7,           blankPlot,
               plot.epsilon.1,       plot.epsilon.2,       plot.epsilon.3,         plot.epsilon.4,
               plot.pi.m,            plot.pi.f,            plot.pi.msm,            blankPlot,
               bottom=textGrob("Priors/posteriors (2/3)", x=0.01, y=1, just="left"),
               ncol=4)

  grid.arrange(plot.rp.1.1.1.1, plot.rp.2.1.1.1, plot.c.min.m1,   plot.c.min.msm1,
               plot.rp.1.1.1.2, plot.rp.2.1.1.2, plot.c.min.m2,   plot.c.min.msm2,
               plot.rp.1.2.1.1, plot.rp.2.2.1.1, plot.c.min.f1,   blankPlot,
               plot.rp.1.2.1.2, plot.rp.2.2.1.2, plot.c.min.f2,   blankPlot,
               plot.rp.1.1.2.1, plot.rp.2.1.2.1, plot.rp.3.1.2.1, plot.rp.4.1.2.1,
               plot.rp.1.1.2.2, plot.rp.2.1.2.2, plot.rp.3.1.2.2, plot.rp.4.1.2.2,
               plot.rp.1.2.2.1, plot.rp.2.2.2.1, plot.rp.3.2.2.1, blankPlot,
               plot.rp.1.2.2.2, plot.rp.2.2.2.2, plot.rp.3.2.2.2, blankPlot,
               bottom=textGrob("Priors/posteriors (3/3)", x=0.01, y=1, just="left"),
               ncol=4)

  dev.off()

}



plot_outputs_from_calibration_file <- function(
  calibration_rds_file,
  site,
  burnin,
  thin = 20,
  sample_size,
  output_directory,
  output_file) {

  if (missing(calibration_rds_file)) stop("The calibration_rds_file path must be provided.")
  if (missing(site)) stop("site must be specified.")
  calibration <- readRDS(calibration_rds_file)
  trace <- as.data.frame(calibration$trace)

  plot_outputs_from_trace(
    site = site,
    trace = trace,
    thin = thin,
    burnin = burnin,
    sample_size = sample_size,
    output_directory = output_directory,
    output_file = output_file)
}


plot_outputs_from_trace <- function(
  site,
  trace,
  thin,
  burnin,
  sample_size,
  output_directory,
  output_file) {

  load_start(site)
  my_trace<-coda::mcmc(trace)
  if (missing(burnin)) burnin <- min(nrow(my_trace) / 2, 50000)
  trace.burn <- fitR::burnAndThin(my_trace, burn=burnin)
  trace.burn.thin <- fitR::burnAndThin(my_trace, thin=thin)
  gc_assign(trace.burn.thin)

  if (missing(sample_size)) sample_size <- min(nrow(trace.burn)/5, 200)
  post.sample <- model_fits(trace.burn, sample.size=sample_size)
  gc_assign(post.sample)
  pred <- as.data.frame(post.sample$outputs)
  gc_assign(pred)
  theta.list <- as.data.frame(post.sample$theta.list)
  if (missing(output_file)) output_file <- paste0(format.Date(Sys.Date(), "%Y-%m-%d"), " ", site, "-calibration.pdf")
  plot_posteriors(output_directory = output_directory, filename = output_file)
}
