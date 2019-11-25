#' Construct a Uniform Likelihood Function
#'
#' In order to come up with trust-able results, we want to
#' try creating a version of the likelihood function
#' which is either 1 in a range around each data point but
#' 0 everywhere outside each range.
#'
#' First, a list of the likelihood components:
#' c(
#'   'Age Assortative Mixing',
#'   'Proportion of Male Cases in MSM 15-24 y',
#'   'Proportion of Male Cases in MSM 25-39 y',
#'   'Proportion of Diagnosed Cases that are Symptomatic',
#'   'Reported Cases',
#'   'Reported Cases Relative Risk',
#'
#'
#' )
#'
#' @return a value 0 or 1
uniform_likelihood <- function(
  pred,
  N = 5
  ) {

  p.symp.ssun <- gc_env$p.symp.ssun
  diag.rate <- gc_env$diag.rate
  rr.diag.subpop <- gc_env$rr.diag.subpop
  rr.diag.subpop.sd <- gc_env$rr.diag.subpop.sd
  p.msm.ssun <- gc_env$p.msm.ssun
  age.dist.dat <- gc_env$age.dist.dat
  var.symp.ssun <- gc_env$var.symp.ssun
  var.p.msm.ssun <- gc_env$var.p.msm.ssun
  diag.rate.sd <- gc_env$diag.rate.sd

  likelihood <- list()

  max_p_symp_sd <- .6/N
  likelihood[['p.symp']] <-
  if(site == "SF") {
    dunif(
      x = as.numeric(unlist(pred[["prev"]][["fit.symp"]]))[2:length(p.symp.ssun)],
      min = (p.symp.ssun - max_p_symp_sd*N)[2:length(p.symp.ssun)],
      max = (p.symp.ssun + max_p_symp_sd*N)[2:length(p.symp.ssun)]
    )
  } else {
    dunif(
      x = as.numeric(unlist(pred[["prev"]][["fit.symp"]])),
      min = p.symp.ssun - max_p_symp_sd*N,
      max = p.symp.ssun + max_p_symp_sd*N
    )
  }

  max_diag_rate_sd <- max(diag.rate.sd)*1.2
  likelihood[['diag.m']] <-
  dunif(
    x = as.numeric(unlist(pred[["prev"]][["fit.diag.rate"]]))[1:(length(diag.rate)/2)],
    min = (diag.rate - max_diag_rate_sd*N)[1:(length(diag.rate)/2)],
    max = (diag.rate + max_diag_rate_sd*N)[1:(length(diag.rate)/2)]
  )

  likelihood[['diag.f']] <-
    dunif(
      x = as.numeric(unlist(pred[["prev"]][["fit.diag.rate"]]))[((length(diag.rate)/2)+1):length(diag.rate)],
      min = (diag.rate - max_diag_rate_sd*N)[((length(diag.rate)/2)+1):length(diag.rate)],
      max = (diag.rate + max_diag_rate_sd*N)[((length(diag.rate)/2)+1):length(diag.rate)]
    )

  max_rr_diag_subpop_sd <- max(rr.diag.subpop.sd)*1.2
  likelihood[['diag.rr']] <-
    dunif(
      x = as.numeric(unlist(pred[["prev"]][["diag.rr"]])),
      min = rr.diag.subpop - max_rr_diag_subpop_sd*N,
      max = rr.diag.subpop + max_rr_diag_subpop_sd*N
    )

  # max_p_msm <- max(sqrt(var.p.msm.ssun))
  max_p_msm <- .25/N
  likelihood[['p.msm']] <-
    dunif(
      x = as.numeric(unlist(pred[["prev"]][["p.diag.msm"]])),
      min = p.msm.ssun - max_p_msm*N,
      max = p.msm.ssun + max_p_msm*N
    )

  max_age_sd <- max(sqrt(age.dist.dat$var))
  likelihood[['age']] <-
    dunif(
      x = pred[["age.dist.all"]][1:4],
      min = (age.dist.dat$p.same.age - max_age_sd*N)[1:4],
      max = (age.dist.dat$p.same.age + max_age_sd*N)[1:4]
    )

  likelihood <- lapply(likelihood, function(x) as.numeric(x != 0))
  return(likelihood)
}

#' Plot Uniform Likelihood Region
#'
#' @return A list of ggplots which represent the likelihood
#' region around each data point.
plot_likelihood_components <- function() {
  require(ggplot2)
  if (! exists('site', envir = gc_env)) stop("site doesn't exist in gc_env")
  if (! exists('plot_environment')) stop("plot_environment doesn't exist")

  with(plot_environment, {
    ### plot diagnosed cases by age and sex ###
    max_diag_rate_sd <- max(diag.rate.sd)*5*1.2
    plot.diag.y.m <- ggplot(data=diag.data[1:length.diag,], aes(x=year, y=diag.rate*100, ymax = (diag.rate + max_diag_rate_sd)*100, ymin = (diag.rate - max_diag_rate_sd)*100))+
      theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
      geom_pointrange(color="red", shape=15, size=.3) +
      labs(title="All M: 15-24y", x="Year", y="Reported cases per 100") +
      scale_x_continuous(breaks = seq((end.year-length.diag+1),end.year,1))

    plot.diag.o.m <- ggplot(data=diag.data[(length.diag+1):(2*length.diag),], aes(x=year, y=diag.rate*100, ymax = (diag.rate + max_diag_rate_sd)*100, ymin = (diag.rate - max_diag_rate_sd)*100))+
      theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
      geom_pointrange(color="red", shape=15, size=.3) +
      labs(title="All M: 25-39y", x="Year", y="Reported cases per 100") +
      scale_x_continuous(breaks = seq((end.year-length.diag+1),end.year,1))

    plot.diag.y.f <- ggplot(data=diag.data[(2*length.diag+1):(3*length.diag),], aes(x=year, y=diag.rate*100, ymax = (diag.rate + max_diag_rate_sd)*100, ymin = (diag.rate - max_diag_rate_sd)*100))+
      theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
      geom_pointrange(color="red", shape=15, size=.3) +
      labs(title="All F: 15-24y", x="Year", y="Reported cases per 100") +
      scale_x_continuous(breaks = seq((end.year-length.diag+1),end.year,1))

    plot.diag.o.f <- ggplot(data=diag.data[(3*length.diag+1):(4*length.diag),], aes(x=year, y=diag.rate*100, ymax = (diag.rate + max_diag_rate_sd)*100, ymin = (diag.rate - max_diag_rate_sd)*100))+
      theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
      geom_pointrange(color="red", shape=15, size=.3) +
      labs(title="All F: 25-39y", x="Year", y="Reported cases per 100") +
      scale_x_continuous(breaks = seq((end.year-length.diag+1),end.year,1))

    max_rr_diag_subpop_sd <- max(rr.diag.subpop.sd)*5
    plot.rr.diag.m.y.1 <- ggplot(data=diag.rr[1:length.rr,], aes(x=year, y=rr_diag, ymin=rr_diag-max_rr_diag_subpop_sd, ymax=rr_diag+max_rr_diag_subpop_sd))+
      geom_pointrange(color="red", shape=15, size=0.3) +
      theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
      labs(title="Black vs all M: 15-24", x="Year", y="Reported case relative risk") +
      scale_x_continuous(breaks = seq((end.year-length.diag+1),end.year,1))

    plot.rr.diag.m.o.1 <- ggplot(data=diag.rr[(length.rr+1):(2*length.rr),],aes(x=year, y=rr_diag, ymin=rr_diag-max_rr_diag_subpop_sd, ymax=rr_diag+max_rr_diag_subpop_sd))+
      geom_pointrange(color="red", shape=15, size=0.3) +
      theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
      labs(title="Black vs all M: 25-39", x="Year", y="Reported case relative risk") +
      scale_x_continuous(breaks = seq((end.year-length.diag+1),end.year,1))

    plot.rr.diag.m.y.2 <- ggplot(data=diag.rr[(2*length.rr+1):(3*length.rr),],aes(x=year, y=rr_diag, ymin=rr_diag-max_rr_diag_subpop_sd, ymax=rr_diag+max_rr_diag_subpop_sd))+
      geom_pointrange(color="red", shape=15, size=0.3) +
      theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
      labs(title="Hispanic vs all M: 15-24", x="Year", y="Reported case relative risk") +
      scale_x_continuous(breaks=seq((end.year-length.rr+1),end.year,1))

    plot.rr.diag.m.o.2 <- ggplot(data=diag.rr[(3*length.rr+1):(4*length.rr),],aes(x=year, y=rr_diag, ymin=rr_diag-max_rr_diag_subpop_sd, ymax=rr_diag+max_rr_diag_subpop_sd))+
      geom_pointrange(color="red", shape=15, size=0.3) +
      theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
      labs(title="Hispanic vs all M: 25-39", x="Year", y="Reported case relative risk") +
      scale_x_continuous(breaks=seq((end.year-length.rr+1),end.year,1))

    plot.rr.diag.f.y.1 <- ggplot(data=diag.rr[(4*length.rr+1):(5*length.rr),],aes(x=year, y=rr_diag, ymin=rr_diag-max_rr_diag_subpop_sd, ymax=rr_diag+max_rr_diag_subpop_sd))+
      geom_pointrange(color="red", shape=15, size=0.3) +
      theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
      labs(title="Black vs all F: 15-24", x="Year", y="Reported case relative risk") +
      scale_x_continuous(breaks=seq((end.year-length.rr+1),end.year,1))

    plot.rr.diag.f.o.1 <- ggplot(data=diag.rr[(5*length.rr+1):(6*length.rr),],aes(x=year, y=rr_diag, ymin=rr_diag-max_rr_diag_subpop_sd, ymax=rr_diag+max_rr_diag_subpop_sd))+
      geom_pointrange(color="red", shape=15, size=0.3) +
      theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
      labs(title="Black vs all F: 25-39", x="Year", y="Reported case relative risk") +
      scale_x_continuous(breaks=seq((end.year-length.rr+1),end.year,1))

    plot.rr.diag.f.y.2 <- ggplot(data=diag.rr[(6*length.rr+1):(7*length.rr),],aes(x=year, y=rr_diag, ymin=rr_diag-max_rr_diag_subpop_sd, ymax=rr_diag+max_rr_diag_subpop_sd))+
      geom_pointrange(color="red", shape=15, size=0.3) +
      theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
      labs(title="Hispanic vs all F: 15-24", x="Year", y="Reported case relative risk") +
      scale_x_continuous(breaks=seq((end.year-length.rr+1),end.year,1))

    plot.rr.diag.f.o.2 <- ggplot(data=diag.rr[(7*length.rr+1):(8*length.rr),],aes(x=year, y=rr_diag, ymin=rr_diag-max_rr_diag_subpop_sd, ymax=rr_diag+max_rr_diag_subpop_sd))+
      geom_pointrange(color="red", shape=15, size=0.3) +
      theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
      labs(title="Hispanic vs all F: 25-39", x="Year", y="Reported case relative risk") +
      scale_x_discrete(labels=seq((end.year-length.rr+1),end.year,1))

    age.dist.dat$label <- c(
      "M young",
      "M old",
      "F young",
      "F old",
      "M young black",
      "M old black",
      "M young other",
      "M old other",
      "M young hispanic",
      "M old hispanic",
      "F young black",
      "F old black",
      "F young other",
      "F old other",
      "F young hispanic",
      "F old hispanic",
      "MSM young",
      "MSM old"
    )

    max_age_sd <- max(sqrt(age.dist.dat$var))*5

    plot.age <- ggplot(data=cbind.data.frame(age.dist.dat, sd = sqrt(age.dist.dat$var)),aes(x=label, y=p.same.age, ymin=p.same.age-max_age_sd, ymax=p.same.age+max_age_sd))+
      geom_pointrange(color="red", shape=15, size=0.3) +
      theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1),title=element_text(size=8)) +
      labs(title="Age assortative mixing", x="Population", y="Proportion of partnerships \nwith same age group")

    max_p_symp_sd <- .6

    plot.p.symp <- ggplot() +
      geom_pointrange(
        data = cbind.data.frame(
          population = c(
            "MSW 15-24 y",
            "MSW 25-39 y",
            "MSM 15-24 y",
            "MSM 25-39 y",
            "F 15-24 y",
            "F 25-39 y"
          ),
          p.symp.dat,
          sd = sqrt(p.symp.dat$p.symp.var)
        ),
        mapping = aes(x = population, y = p.symp, ymin = p.symp - max_p_symp_sd, ymax = p.symp + max_p_symp_sd),
        color = "red",
        shape = 15,
        size = .3
      ) +
      theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
      labs(title="Proportion of diagnosed cases\nthat are symptomatic", x="Population", y="Proportion of cases that are symptomatic")

    max_p_msm <- .5
    plot.p.msm.y <-
      ggplot(
        data = as.data.frame(p.msm.plot[1:length.msm]),
        mapping = aes(
          x = seq(end.year - length.msm + 1, end.year, 1),
          y = p.msm.plot[1:length.msm],
          ymax = p.msm.plot[1:length.msm] + max_p_msm,
          ymin = p.msm.plot[1:length.msm] - max_p_msm
        )) +
      geom_pointrange(color = "red",
                 shape = 15,
                 size = .3) +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1),
        title = element_text(size = 8)
      ) +
      labs(title = "Proportion of male cases\n diagnosed in MSM 15-24 y", x =
             "Year", y = "Proportion of cases in MSM") +
      scale_x_continuous(breaks = seq(end.year - length.msm + 1, end.year, 1))

    plot.p.msm.o <-
      ggplot(
        data = as.data.frame(p.msm.plot[(length.msm + 1):(2 * length.msm)]),
        aes(x = seq((end.year - length.msm + 1), end.year, 1),
            y = p.msm.plot[(length.msm + 1):(2 * length.msm)],
            ymax = p.msm.plot[(length.msm + 1):(2 * length.msm)] + max_p_msm,
            ymin = p.msm.plot[(length.msm + 1):(2 * length.msm)] - max_p_msm
        )) +
      geom_pointrange(color = "red",
                 shape = 15,
                 size = .3) +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1),
        title = element_text(size = 8)
      ) +
      labs(title = "Proportion of male cases\n diagnosed in MSM 25-39y y", x =
             "Year", y = "Proportion of cases in MSM") +
      scale_x_continuous(breaks = seq((end.year - length.msm + 1), end.year, 1))

  })
  return(invisible(NULL))
}

#' Prepare plot_environment with necessary computations and values
prepare_for_plotting_likelihoods <- function() {
  require(reshape2)
  plot_environment <- new.env()
  with(plot_environment, {
    # attach(gc_env)
    for(n in ls(gc_env, all.names=TRUE)) assign(n, get(n, gc_env), plot_environment)
    # prep input data
    diag.age.sex.rate <- if(!"year" %in% colnames(as.data.frame(diag.age.sex.rate))) tibble::rownames_to_column(as.data.frame(diag.age.sex.rate), "year")  else diag.age.sex.rate
    diag.dat <- melt(diag.age.sex.rate[,1:5],  id.vars = "year", measured.vars=c("m_y", "m_o", "f_y","f_o"))
    diag.dat$cat <- ifelse(diag.dat$variable=="m_y", "M 15-24 y", ifelse(diag.dat$variable=="m_o", "M 25-39 y", ifelse(diag.dat$variable=="f_y", "F 15-24 y", "F 25-39 y")))
    p.msm.dat$cat <- ifelse(p.msm.dat$age_cat==1, "15-24 y", "25-39 y" )
    diag.subpop.rate <- if(!"year" %in% colnames(diag.subpop.rate)) tibble::rownames_to_column(diag.subpop.rate, "year")  else diag.subpop.rate
    subpop.dat <- melt(diag.subpop.rate,  id.vars = "year", measured.vars=list[2:ncol(diag.subpop.rate)])
    subpop.dat$cat <- ifelse(
      subpop.dat$variable == "y.m.1", "Black M 15-24 y",
      ifelse(
        subpop.dat$variable == "y.m.2",
        "Other M 15-24 y",
        ifelse(
          subpop.dat$variable == "y.m.3",
          "Hispanic M 15-24 y",
          ifelse(
            subpop.dat$variable == "o.m.1",
            "Black M 25-39 y",
            ifelse(
              subpop.dat$variable == "o.m.2",
              "Other M 25-39 y",
              ifelse(
                subpop.dat$variable == "o.m.3",
                "Hispanic M 25-39 y",
                ifelse(
                  subpop.dat$variable == "y.f.1",
                  "Black F 15-24 y",
                  ifelse(
                    subpop.dat$variable == "y.f.2",
                    "Other F 15-24 y",
                    ifelse(
                      subpop.dat$variable == "y.f.3",
                      "Hispanic F 15-24 y",
                      ifelse(
                        subpop.dat$variable == "o.f.1",
                        "Black F 25-39 y",
                        ifelse(
                          subpop.dat$variable == "o.f.2",
                          "Other F 25-39 y",
                          "Hispanic F 25-39 y"
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      ))

    subpop.dat$sex <- ifelse(subpop.dat$variable=="y.m.1"|subpop.dat$variable=="y.m.2"|subpop.dat$variable=="y.m.3"|
                               subpop.dat$variable=="o.m.1"|subpop.dat$variable=="o.m.2"|subpop.dat$variable=="o.m.3", "M", "F")


    subpop.dat$age <- ifelse(subpop.dat$variable=="y.m.1"|subpop.dat$variable=="y.m.2"|subpop.dat$variable=="y.m.3"|
                               subpop.dat$variable=="y.f.1"|subpop.dat$variable=="y.f.2"|subpop.dat$variable=="y.f.3", "15-24 y", "25-39 y")

    subpop.dat$Population <- ifelse(subpop.dat$variable=="y.m.1"|subpop.dat$variable=="o.m.1"|subpop.dat$variable=="y.f.1"|subpop.dat$variable=="o.f.1", "Black",
                                    ifelse(subpop.dat$variable=="y.m.2"|subpop.dat$variable=="o.m.2"|subpop.dat$variable=="y.f.2"|subpop.dat$variable=="o.f.2", "Other", "Hispanic"))

    length.rr <- nrow(diag.rr)/8
    length.msm <- nrow(p.msm.dat)/2
    length.diag <- nrow(diag.age.sex.rate)

    diag.data <- cbind(year=(seq(as.numeric(diag.age.sex.rate[1,"year"]),as.numeric(diag.age.sex.rate[nrow(diag.age.sex.rate), "year"]),1)), as.data.frame(diag.rate), cat=rep(c("y.m", "o.m", "y.f", "o.f"),each=nrow(diag.age.sex.rate) ))
    diag.data$sd <- diag.rate.sd
    age.data <- as.data.frame(age.dist.dat[1:4,])
    p.msm.plot <- p.msm.ssun

    max.y.m <- 100* max(diag.subpop.rate[,"y.m.1"],diag.subpop.rate[,"y.m.2"], diag.subpop.rate[,"y.m.3"], na.rm=TRUE)
    max.o.m <- 100*max(diag.subpop.rate[,"o.m.2"],diag.subpop.rate[,"o.m.3"], na.rm=TRUE)
    max.y.f <- 100*max(diag.subpop.rate[,"y.f.1"], diag.subpop.rate[,"y.f.2"],diag.subpop.rate[,"y.f.3"], na.rm=TRUE)
    max.o.f <- 100*max(diag.subpop.rate[,"o.f.1"],diag.subpop.rate[,"o.f.2"], diag.subpop.rate[,"o.f.3"], na.rm=TRUE)
    max.diag <- max(max.y.m, max.o.m, max.y.f, max.o.f) +0.1
    # max.msm <- max(out.diag.y.msm$value,out.diag.o.msm$value, na.rm=TRUE)*100+0.1
    max.rr.diag <- max(diag.rr$rr_diag*1.2, na.rm=TRUE)
    max.rr.diag.1 <- max(diag.rr$rr_diag*1.2, na.rm=TRUE)+1
    max.rr.diag.2 <- max(diag.rr$rr_diag*1.2, na.rm=TRUE)+0.5



  })
  return(plot_environment)
}


#' Save Plots from plot_environment
save_plots_from_plot_environment <- function() {
  plot_list <- ls(name = plot_environment)
  plot_list <- grep(plot_list, pattern = "^plot", value = T)

  for (p in plot_list) {
    plt <- get(p, envir = plot_environment)
    ggsave(plt, filename = paste0("~/Desktop/Likelihoods/", gsub("\\.", "_", p), ".png"))
  }
}


sample_priors_for_uniform_likelihood_acceptable_thetas <- function(N, out_dir) {
  for (j in 1:N) {
    theta_list <- list()
    for (i in 1:1000) {
      try({
        theta <- sample_prior()
        pred <- prediction_epi(theta)
        out <- uniform_likelihood(pred, 15)
        success <- all(sapply(out, function(x) all(as.logical(x))))
        if (success) {
          theta_list[[length(theta_list)+1]] <- theta
        }
      })
    }
    saveRDS(do.call(rbind, theta_list), file.path(out_dir, paste0("theta_list_", gc_env$site, "_", j, ".rds")))
  }
}
