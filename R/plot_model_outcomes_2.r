###########################################################################
### visualize prior and posterior distributions and calibration targets ###
###########################################################################

# This is a 2nd draft version of the plot_posteriors function file.
# We will add sources to the data targets plots and legends detailing the
# data depicted.



#' Plot Calibration Targets
#'
#' Construct a named list of ggplot objects, each of which depicts one of the
#' pairs of data targets and the model fit against it.
#'
#' @export
#' @import ggplot2
#' @import cowplot
list_target_and_fit_plots <- function() {
  if (! exists('pred')) stop("The pred data.frame does not exist.")

  # reshape the original data and the simulations data in gc_env and pred to
  # work with ggplot2
  prepare_data_for_prior_posterior_target_and_fit_plots()

  # construct the NHANES estimates and comparison
  with(gc_env, {

    { # NHANES Targets ----
    plot.prev <- ggplot(data = out.prev) +
      # data targets
			geom_pointrange(data = filter(nhanes.updated.dat, sex=='W') %>%
											merge(x = ., y = data.frame(shortname = names(longname), longname =
																			 longname), by='shortname'),
                      aes(
                        x = longname,
                        y = prev,
                        ymin = prev - prev_std,
                        ymax = prev + prev_std,
                        color = 'red'
                      ), shape = 15, alpha = 0.7) +
      # model simulations summarized
      geom_boxplot(
        aes(
          x = variable,
          y = value,
          fill = variable,
          alpha = 'boxes'
        ),
        width = 0.75,
        outlier.alpha = 0.3,
        position = "dodge"
      ) +
      # blank theme
      theme_classic() +
      # put the legend centered-right
      theme(
        legend.position=c(.78,.5),
        legend.direction = 'vertical',
        legend.spacing.y = unit(0.1, 'in'),
        # legend.position = 'right',
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)) +
        # legend.position = 'bottom') +
      # use color and alpha as our legend features for NHANES vs. Model Fits
      scale_color_manual(
        values = c('red' = 'red', 'std_clinic' = 'cornflowerblue'),
        label = c("NHANES Estimates", "Clinic Estimate\n(Not Calibrated To)")) +
      scale_alpha_manual(values = 0.6, label="By Population") +
      # labels, re-label color and alpha for our legend
      labs(
        title = "Gonorrhea Prevalence",
        x = "Population",
        y = "Prevalence (%)",
        color = 'Data Targets',
        alpha = 'Model Fits',
        shape = 'STD Clinic Estimate'
      ) +
      # flip coordinates for readability of demographic groups
      coord_flip() +
      # order the x-axis by the levels assigned in out.prev inside
      # prepare_data_for_prior_posterior_target_and_fit_plots
      scale_x_discrete(limits = rev(levels(out.prev$variable))) +
      # expand the quantitative axis by x1.05 so that text isn't clipped off
      expand_limits(y = c(0, max(out.prev$value, nhanes.updated.dat$prev + nhanes.updated.dat$prev_std)*1.05)) +
      # title
      ggtitle("Gonorrhea Prevalence") +
      # turn off fill = demographic group legend.
      guides(fill = F)

      # If the site is San Francisco, include the STD clinic estimate
      if (gc_env$site == "SF") {
        plot.prev <- plot.prev +
          geom_pointrange(
            data = data.frame(x = "All MSM", y = .06),
            aes(
              x = x,
              y = y,
              ymax = y + 0.012,
              ymin = y - 0.012,
              color = 'std_clinic'
            ),
            # color = "dodgerblue",
            shape = 15,
            alpha = .7
          )
      }
    }

    { # Subpopulation assortative mixing
      plot.subpop <- ggplot(data=out.subpop)+
        geom_boxplot(aes(x=Var2, y=value, fill=Var2, alpha = 'subpop'), outlier.alpha = 0.3) +
        geom_pointrange(data=s.dist.sd,aes(x=1:nrow(s.dist.sd), y=dat.s.dist,ymin=lcl,ymax=ucl, color = 'data'), shape=15, size=0.5, alpha=0.8) +
        theme_classic() +
        theme(
          # legend.position = "none",
          title = element_text(size = 8),
          axis.text.x = element_text(angle = 90),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7)
        ) +
        labs(
          title = "Subpopulation assortative mixing\n(not used in calibration)",
          x = "",
          y = "Proportion of partnerships \nwith same subpopulation",
          alpha = 'Model Outcomes',
          color = 'Comparison Data') +
        scale_x_discrete(
          labels = c(
            "Male Black",
            "Male Other",
            "Male Hispanic",
            "Female Black",
            "Female Other",
            "Female Hispanic")) +
        scale_alpha_manual(values = 1, label = '') +
        scale_color_manual(values = 'dimgrey', label = 'NSFG') +
        guides(fill = F) +
        expand_limits(y = c(0,1))
    }

    { # Age Assortativity
      plot.age <- ggplot(data=out.age)+
        geom_boxplot(aes(x=Var2, y=value, fill=Var2, alpha = 'age'), outlier.alpha = 0.3) +
        geom_pointrange(data=age.data,aes(x=1:nrow(age.data), y=p.same.age, ymin=lcl, ymax=ucl, color = 'data'), shape=15, size=0.5) +
        theme_classic() +
        theme(
          title = element_text(size = 8),
          axis.text.x = element_text(angle = 90),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7)
        ) +
        labs(title = "Age assortative mixing",
             x = "Population",
             y = "Proportion of partnerships \nwith same age group",
             alpha = 'Model Fits',
             color = 'Data Targets') +
        scale_alpha_manual(values = 1, label = '') +
        scale_color_manual(values = 'red', label = 'NSFG') +
        scale_x_discrete(labels=c("Male 15-24y","Male 25-39y","Female 15-24y", "Female 25-39y")) +
        guides(fill = F) +
        expand_limits(y = c(0,1))
        # coord_flip(ylim = c(0,1))
    }

    { # Proportion Symptomatic

      color_scale_psymp <- if (gc_env$site == 'SF') c('red',rep('dimgrey', length=nrow(out.symp)-1)) else rep('red', length=nrow(out.symp))
      shape_scale_psymp <- if (gc_env$site == 'SF') c(15, rep(4, length = nrow(out.symp)-1)) else rep(15, length=nrow(out.symp))

      data_labels_psymp <- if (gc_env$site == 'SF') c('NSFG', 'Not Calibrated To\nDue To Only One Case') else 'NSFG'

      plot.p.symp <- ggplot(data=out.symp)+
        geom_boxplot(aes(x=reorder(Var2, desc(Var2)), y=value, fill=Var2, alpha = 'fits'), outlier.alpha = 0.3) +
        geom_pointrange(data=p.symp.dat,
                        aes(x=rev(1:length(p.symp.ssun)), y=p.symp,
                            ymax = p.symp + sqrt(p.symp.var),
                            ymin = p.symp - sqrt(p.symp.var),
                            color = if (gc_env$site == 'SF') as.factor(ifelse(sex == 'M' & age_cat==1, 1, 0)) else '1',
                            shape = if (gc_env$site == 'SF') as.factor(ifelse(sex == 'M' & age_cat==1, 1, 0)) else '1'),
                        size=.5) +
        theme_classic() +
        theme(
          # axis.title.y = element_text(margin = margin(t = 0, r = -60, b = 0, l = 0)),
          title = element_text(size = 8),
          legend.position = 'bottom',
          # legend.direction = 'vertical',
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7)
        ) +
        labs(title = "Proportion of diagnoses\nthat are symptomatic", x =
               "\nPopulation", y = "Proportion of diagnoses that are symptomatic",
             shape = "Data Targets", color = 'Data Targets', alpha = 'Model Fits') +
        scale_x_discrete(# limits = rev(levels(out.symp$Var2)),
          labels=rev(setNames(c("MSW 15-24 y","MSW 25-39 y", "MSM 15-24 y", "MSM 25-39 y","Female 15-24 y", "Female 25-39 y"), levels(out.symp$Var2) ))) +
        scale_color_manual(values=color_scale_psymp, labels = data_labels_psymp) +
        scale_shape_manual(values=shape_scale_psymp, labels = data_labels_psymp) +
        scale_alpha_manual(values = 1, label = "") +
        guides(fill = F) +
        coord_flip(ylim = c(0,1))

    }

    { # Proportion of Cases in MSM
      plot.p.msm.y <- ggplot(data = out.p.msm.y)+
        geom_line(aes(x=Var2, y = value, group = Var1, color = 'grey')) +
        stat_summary(fun.y = mean, geom='line', aes(x = Var2, y=value, group=1, color = 'black')) +
        geom_pointrange(
          data = as.data.frame(p.msm.plot[1:length.msm]),
          aes(x = 1:length.msm,
              y = p.msm.plot[1:length.msm],
              ymax = p.msm.plot[1:length.msm] + sqrt(p.msm.dat$var[1:length.msm]),
              ymin = p.msm.plot[1:length.msm] - sqrt(p.msm.dat$var[1:length.msm]),
              alpha = 'data'),
          color = "red",
          shape = 15,
          size = .5
        ) +
        theme_classic() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1),
          title = element_text(size = 8),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7)
        ) +
        labs(title = "Proportion of male diagnoses\n in MSM 15-24y", x =
               "Year", y = "Proportion of diagnoses in MSM", alpha = 'Data Targets', color = 'Model Fits') +
        coord_cartesian(ylim=c(0,1))+
        scale_x_discrete(labels=seq(end.year-length.msm + 1, end.year, 1)) +
        scale_color_manual(values = c('grey' = 'grey', 'black' = 'black'), labels = c('grey' = 'Individual Simulations', 'black' = 'Mean of Simulations')) +
        scale_alpha_manual(values = 1, labels = 'SSuN')

      plot.p.msm.o <- ggplot(data=out.p.msm.o)+
        geom_line(aes(x=Var2, y = value, group = Var1, color = 'grey')) +
        stat_summary(fun.y = mean, geom='line', aes(x = Var2, y=value, group=1, color = 'black')) +
        geom_pointrange(
          data = as.data.frame(p.msm.plot[(length.msm + 1):(2 * length.msm)]),
          aes(
            x = 1:length.msm,
            y = p.msm.plot[(length.msm + 1):(2 * length.msm)],
            ymax = p.msm.plot[(length.msm + 1):(2 * length.msm)] + sqrt(p.msm.dat$var[(length.msm + 1):(2 * length.msm)]),
            ymin = p.msm.plot[(length.msm + 1):(2 * length.msm)] - sqrt(p.msm.dat$var[(length.msm + 1):(2 * length.msm)]),
            alpha = 'data'),
          color = "red",
          shape = 15,
          size = .5
        ) +
        theme_classic() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1),
          title = element_text(size = 8),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7)
        ) +
        labs(title = "Proportion of male diagnoses\n in MSM 25-39y", x =
               "Year", y = "Proportion of diagnoses in MSM", alpha = 'Data Targets', color = 'Model Fits') +
        coord_cartesian(ylim=c(0,1))+
        scale_x_discrete(labels=seq((end.year-length.msm+1),end.year,1)) +
        scale_color_manual(values = c('grey' = 'grey', 'black' = 'black'), labels = c('grey' = 'Individual Simulations', 'black' = 'Mean of Simulations')) +
        scale_alpha_manual(values = 1, labels = 'SSuN')

    }

    gc_env$target_and_fit_plots_list <-
      list(plot.prev = plot.prev,
           plot.subpop = plot.subpop,
           plot.age = plot.age,
           plot.p.symp = plot.p.symp,
           plot.p.msm.y = plot.p.msm.y,
           plot.p.msm.o = plot.p.msm.o
           )
  })


  # return list of plots
  gc_env$target_and_fit_plots_list
}

#' Prepare Data for Prior, Posterior, Target, and Model Fit Plots
#'
#' @import reshape2
prepare_data_for_prior_posterior_target_and_fit_plots <- function() {
  with(gc_env, {

    #prep input data
    diag.age.sex.rate <-
      if (!"year" %in% colnames(as.data.frame(diag.age.sex.rate)))
        tibble::rownames_to_column(as.data.frame(diag.age.sex.rate), "year")
    else
      diag.age.sex.rate
    diag.dat <-
      melt(
        diag.age.sex.rate[, 1:5],
        id.vars = "year",
        measured.vars = c("m_y", "m_o", "f_y", "f_o")
      )
    diag.dat$cat <-
      ifelse(
        diag.dat$variable == "m_y",
        "M 15-24 y",
        ifelse(
          diag.dat$variable == "m_o",
          "M 25-39 y",
          ifelse(diag.dat$variable == "f_y", "F 15-24 y", "F 25-39 y")
        )
      )
    p.msm.dat$cat <-
      ifelse(p.msm.dat$age_cat == 1, "15-24 y", "25-39 y")
    diag.subpop.rate <-
      if (!"year" %in% colnames(diag.subpop.rate))
        tibble::rownames_to_column(diag.subpop.rate, "year")
    else
      diag.subpop.rate
    subpop.dat <-
      melt(diag.subpop.rate,
           id.vars = "year",
           measured.vars = list[2:ncol(diag.subpop.rate)])
    subpop.dat$cat <-
      ifelse(
        subpop.dat$variable == "y.m.1",
        "Black M 15-24 y",
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
        )
      )

    subpop.dat$sex <-
      ifelse(
        subpop.dat$variable == "y.m.1" |
          subpop.dat$variable == "y.m.2" | subpop.dat$variable == "y.m.3" |
          subpop.dat$variable == "o.m.1" |
          subpop.dat$variable == "o.m.2" |
          subpop.dat$variable == "o.m.3",
        "M",
        "F"
      )


    subpop.dat$age <-
      ifelse(
        subpop.dat$variable == "y.m.1" |
          subpop.dat$variable == "y.m.2" | subpop.dat$variable == "y.m.3" |
          subpop.dat$variable == "y.f.1" |
          subpop.dat$variable == "y.f.2" |
          subpop.dat$variable == "y.f.3",
        "15-24 y",
        "25-39 y"
      )

    subpop.dat$Population <-
      ifelse(
        subpop.dat$variable == "y.m.1" |
          subpop.dat$variable == "o.m.1" |
          subpop.dat$variable == "y.f.1" |
          subpop.dat$variable == "o.f.1",
        "Black",
        ifelse(
          subpop.dat$variable == "y.m.2" |
            subpop.dat$variable == "o.m.2" |
            subpop.dat$variable == "y.f.2" |
            subpop.dat$variable == "o.f.2",
          "Other",
          "Hispanic"
        )
      )

    # prep outputs for plotting
    # here y=15-24y, o=25-39y, m=male, f=female, 1=black, 2=other, 3=Hispanic, msm=men who have sex with men
    length.rr <- nrow(diag.rr) / 8
    # out.prev <- melt(as.matrix(subset(pred, select=prev.fit.prev.extra1:prev.fit.prev.extra15)))  #model prevalence with extra categories, with MSM in calculations
    out.prev <-
      pred[, grep(
        'prev.msm_redistributed_prevalence_rates|prev.fit.prev.extra.prev.msm',
        names(pred),
        value = T
      )]
    out.prev <- out.prev[complete.cases(out.prev), ]
    out.diag.y.m.1 <-
      melt(as.matrix(subset(
        pred, select = prev.fit.diag.subpop1:(prev.fit.diag.subpop1 + length.rr -
                                                1)
      ))) #reported cases by subpop
    out.diag.o.m.1 <-
      melt(as.matrix(subset(
        pred,
        select = (prev.fit.diag.subpop1 + length.rr):(prev.fit.diag.subpop1 + 2 *
                                                        length.rr - 1)
      )))
    out.diag.y.m.2 <-
      melt(as.matrix(subset(
        pred,
        select = (prev.fit.diag.subpop1 + 2 * length.rr):(prev.fit.diag.subpop1 +
                                                            3 * length.rr - 1)
      )))
    out.diag.o.m.2 <-
      melt(as.matrix(subset(
        pred,
        select = (prev.fit.diag.subpop1 + 3 * length.rr):(prev.fit.diag.subpop1 +
                                                            4 * length.rr - 1)
      )))
    out.diag.y.m.3 <-
      melt(as.matrix(subset(
        pred,
        select = (prev.fit.diag.subpop1 + 4 * length.rr):(prev.fit.diag.subpop1 +
                                                            5 * length.rr - 1)
      )))
    out.diag.o.m.3 <-
      melt(as.matrix(subset(
        pred,
        select = (prev.fit.diag.subpop1 + 5 * length.rr):(prev.fit.diag.subpop1 +
                                                            6 * length.rr - 1)
      )))
    out.diag.y.f.1 <-
      melt(as.matrix(subset(
        pred,
        select = (prev.fit.diag.subpop1 + 6 * length.rr):(prev.fit.diag.subpop1 +
                                                            7 * length.rr - 1)
      )))
    out.diag.o.f.1 <-
      melt(as.matrix(subset(
        pred,
        select = (prev.fit.diag.subpop1 + 7 * length.rr):(prev.fit.diag.subpop1 +
                                                            8 * length.rr - 1)
      )))
    out.diag.y.f.2 <-
      melt(as.matrix(subset(
        pred,
        select = (prev.fit.diag.subpop1 + 8 * length.rr):(prev.fit.diag.subpop1 +
                                                            9 * length.rr - 1)
      )))
    out.diag.o.f.2 <-
      melt(as.matrix(subset(
        pred,
        select = (prev.fit.diag.subpop1 + 9 * length.rr):(prev.fit.diag.subpop1 +
                                                            10 * length.rr - 1)
      )))
    out.diag.y.f.3 <-
      melt(as.matrix(subset(
        pred,
        select = (prev.fit.diag.subpop1 + 10 * length.rr):(prev.fit.diag.subpop1 +
                                                             11 * length.rr - 1)
      )))
    out.diag.o.f.3 <-
      melt(as.matrix(subset(
        pred,
        select = (prev.fit.diag.subpop1 + 11 * length.rr):(prev.fit.diag.subpop1 +
                                                             12 * length.rr - 1)
      )))
    out.diag.y.msm <-
      melt(as.matrix(subset(
        pred, select = prev.fit.diag.msm1:(prev.fit.diag.msm1 + length.rr - 1)
      )))
    out.diag.o.msm <-
      melt(as.matrix(subset(
        pred,
        select = (prev.fit.diag.msm1 + length.rr):(prev.fit.diag.msm1 + 2 * length.rr -
                                                     1)
      )))

    length.diag <- nrow(diag.age.sex.rate)
    out.diag.y.m <-
      melt(as.matrix(subset(
        pred, select = prev.fit.diag.rate1:(prev.fit.diag.rate1 + length.diag -
                                              1)
      ))) # reported cases rates by age and sex only
    out.diag.o.m <-
      melt(as.matrix(subset(
        pred,
        select = (prev.fit.diag.rate1 + length.diag):(prev.fit.diag.rate1 + 2 *
                                                        length.diag - 1)
      )))
    out.diag.y.f <-
      melt(as.matrix(subset(
        pred,
        select = (prev.fit.diag.rate1 + 2 * length.diag):(prev.fit.diag.rate1 +
                                                            3 * length.diag - 1)
      )))
    out.diag.o.f <-
      melt(as.matrix(subset(
        pred,
        select = (prev.fit.diag.rate1 + 3 * length.diag):(prev.fit.diag.rate1 +
                                                            4 * length.diag - 1)
      )))

    out.inc.y.m.1 <-
      melt(100 * as.matrix(subset(pred, select = prev.inc1:(
        prev.inc1 + length.diag - 1
      ))) / sum(n.i[y.m.1])) #model incidence by age, sex, and subpopulation
    out.inc.o.m.1 <-
      melt(100 * as.matrix(subset(
        pred,
        select = (prev.inc1 + length.diag):(prev.inc1 + 2 * length.diag - 1)
      )) / sum(n.i[o.m.1]))
    out.inc.y.m.2 <-
      melt(100 * as.matrix(subset(
        pred,
        select = (prev.inc1 + 2 * length.diag):(prev.inc1 + 3 * length.diag - 1)
      )) / sum(n.i[y.m.2]))
    out.inc.o.m.2 <-
      melt(100 * as.matrix(subset(
        pred,
        select = (prev.inc1 + 3 * length.diag):(prev.inc1 + 4 * length.diag - 1)
      )) / sum(n.i[o.m.2]))
    out.inc.y.m.3 <-
      melt(100 * as.matrix(subset(
        pred,
        select = (prev.inc1 + 4 * length.diag):(prev.inc1 + 5 * length.diag - 1)
      )) / sum(n.i[y.m.3]))
    out.inc.o.m.3 <-
      melt(100 * as.matrix(subset(
        pred,
        select = (prev.inc1 + 5 * length.diag):(prev.inc1 + 6 * length.diag - 1)
      )) / sum(n.i[o.m.3]))
    out.inc.y.msm <-
      melt(100 * as.matrix(subset(
        pred,
        select = (prev.inc1 + 6 * length.diag):(prev.inc1 + 7 * length.diag - 1)
      )) / sum(n.i[y.m.4]))
    out.inc.o.msm <-
      melt(100 * as.matrix(subset(
        pred,
        select = (prev.inc1 + 7 * length.diag):(prev.inc1 + 8 * length.diag - 1)
      )) / sum(n.i[o.m.4]))
    out.inc.y.f.1 <-
      melt(100 * as.matrix(subset(
        pred,
        select = (prev.inc1 + 8 * length.diag):(prev.inc1 + 9 * length.diag - 1)
      )) / sum(n.i[y.f.1]))
    out.inc.o.f.1 <-
      melt(100 * as.matrix(subset(
        pred,
        select = (prev.inc1 + 9 * length.diag):(prev.inc1 + 10 * length.diag -
                                                  1)
      )) / sum(n.i[o.f.1]))
    out.inc.y.f.2 <-
      melt(100 * as.matrix(subset(
        pred,
        select = (prev.inc1 + 10 * length.diag):(prev.inc1 + 11 * length.diag -
                                                   1)
      )) / sum(n.i[y.f.2]))
    out.inc.o.f.2 <-
      melt(100 * as.matrix(subset(
        pred,
        select = (prev.inc1 + 11 * length.diag):(prev.inc1 + 12 * length.diag -
                                                   1)
      )) / sum(n.i[o.f.2]))
    out.inc.y.f.3 <-
      melt(100 * as.matrix(subset(
        pred,
        select = (prev.inc1 + 12 * length.diag):(prev.inc1 + 13 * length.diag -
                                                   1)
      )) / sum(n.i[y.f.3]))
    out.inc.o.f.3 <-
      melt(100 * as.matrix(subset(
        pred,
        select = (prev.inc1 + 13 * length.diag):(prev.inc1 + 14 * length.diag -
                                                   1)
      )) / sum(n.i[o.f.3]))

    out.rr.diag.m.y.1 <-
      melt(as.matrix(subset(
        pred, select = prev.diag.rr1:(prev.diag.rr1 + length.rr - 1)
      ))) #relative risk of being a reported case by subpop, age, sex (here 1=black, 2=Hispanic)
    out.rr.diag.m.o.1 <-
      melt(as.matrix(subset(
        pred,
        select = (prev.diag.rr1 + length.rr):(prev.diag.rr1 + 2 * length.rr - 1)
      )))
    out.rr.diag.m.y.2 <-
      melt(as.matrix(subset(
        pred,
        select = (prev.diag.rr1 + 2 * length.rr):(prev.diag.rr1 + 3 * length.rr -
                                                    1)
      )))
    out.rr.diag.m.o.2 <-
      melt(as.matrix(subset(
        pred,
        select = (prev.diag.rr1 + 3 * length.rr):(prev.diag.rr1 + 4 * length.rr -
                                                    1)
      )))
    out.rr.diag.f.y.1 <-
      melt(as.matrix(subset(
        pred,
        select = (prev.diag.rr1 + 4 * length.rr):(prev.diag.rr1 + 5 * length.rr -
                                                    1)
      )))
    out.rr.diag.f.o.1 <-
      melt(as.matrix(subset(
        pred,
        select = (prev.diag.rr1 + 5 * length.rr):(prev.diag.rr1 + 6 * length.rr -
                                                    1)
      )))
    out.rr.diag.f.y.2 <-
      melt(as.matrix(subset(
        pred,
        select = (prev.diag.rr1 + 6 * length.rr):(prev.diag.rr1 + 7 * length.rr -
                                                    1)
      )))
    out.rr.diag.f.o.2 <-
      melt(as.matrix(subset(
        pred,
        select = (prev.diag.rr1 + 7 * length.rr):(prev.diag.rr1 + 8 * length.rr -
                                                    1)
      )))

    out.subpop <-
      melt(as.matrix(subset(pred, select = pred.s.dist1:pred.s.dist6))) #model subpopulation assortativity
    out.age <-
      melt(as.matrix(subset(pred, select = age.dist.all1:age.dist.all4))) #model age assortativity
    out.symp <-
      melt(as.matrix(subset(pred, select = prev.fit.symp1:prev.fit.symp6))) #model proportion of diagnosed cases that are symptomatic
    length.msm <- nrow(p.msm.dat) / 2
    out.p.msm.y <-
      melt(as.matrix(subset(
        pred, select = prev.p.diag.msm1:(prev.p.diag.msm1 + length.msm - 1)
      )))
    out.p.msm.o <-
      melt(as.matrix(subset(
        pred,
        select = (prev.p.diag.msm1 + length.msm):(prev.p.diag.msm1 + 2 * length.msm -
                                                    1)
      )))

    diag.data <-
      cbind(year = (seq(
        as.numeric(diag.age.sex.rate[1, "year"]),
        as.numeric(diag.age.sex.rate[nrow(diag.age.sex.rate), "year"]),
        1
      )),
      as.data.frame(diag.rate),
      cat = rep(c("y.m", "o.m", "y.f", "o.f"), each = nrow(diag.age.sex.rate)))
    age.data <- as.data.frame(age.dist.dat[1:4, ])
    p.msm.plot <- p.msm.ssun

    nhanes.updated.dat <- gc_env$nhanes.updated.dat

    shortname <-
      c(
        "Y F",
        "Y M",
        "Y F Black",
        "Y M Black",
        "Y F Other",
        "Y M Other",
        "Y F Hisp",
        "Y M Hisp",
        "O F",
        "O M",
        "O F Black",
        "O M Black",
        "O F Other",
        "O M Other",
        "O F Hisp",
        "O M Hisp"
      )

    longname <- sapply(shortname, function(x) {
      x <- strsplit(x, " ")[[1]]
      lookup_names <- c(
        Y = "15-24y",
        O = "25-39y",
        F = "Female",
        M = "Male",
        Black = "Black",
        Other = "Other",
        Hisp = "Hispanic",
        All = "All",
        MSM = "MSM"
      )
      paste0(sapply(x, function(x1) lookup_names[x1]), collapse = " ")
    })

    nhanes_cols <-
      grepl("prev.msm_redistributed_prevalence_rates.prev",
            colnames(out.prev))
    colnames(out.prev)[nhanes_cols] <- longname
    colnames(out.prev) <-
      gsub("prev.fit.prev.extra.prev.msm", "All MSM", colnames(out.prev))
    nhanes.updated.dat[['shortname']] <- shortname

    out.prev <- reshape2:::melt.data.frame(out.prev, measure.vars = colnames(out.prev))

    out.prev$variable <- factor(out.prev$variable, ordered = T, levels = c(

      "15-24y Female",
      "15-24y Male",

      "25-39y Female",
      "25-39y Male",

      "15-24y Female Black",
      "15-24y Female Hispanic",
      "15-24y Female Other",

      "15-24y Male Black",
      "15-24y Male Hispanic",
      "15-24y Male Other",

      "25-39y Female Black",
      "25-39y Female Hispanic",
      "25-39y Female Other",

      "25-39y Male Black",
      "25-39y Male Hispanic",
      "25-39y Male Other",

      "All MSM"
    ))



    #get max values of reported cases for setting y-axis
    max.y.m <-
      100 * max(
        out.diag.y.m.1$value,
        out.diag.y.m.2$value,
        out.diag.y.m.3$value,
        diag.subpop.rate[, "y.m.1"],
        diag.subpop.rate[, "y.m.2"],
        diag.subpop.rate[, "y.m.3"],
        na.rm = TRUE
      )
    max.o.m <-
      100 * max(
        out.diag.o.m.1$value,
        out.diag.o.m.2$value,
        out.diag.o.m.3$value,
        diag.subpop.rate[, "o.m.1"],
        diag.subpop.rate[, "o.m.2"],
        diag.subpop.rate[, "o.m.3"],
        na.rm = TRUE
      )
    max.y.f <-
      100 * max(
        out.diag.y.f.1$value,
        out.diag.y.f.2$value,
        out.diag.y.f.3$value,
        diag.subpop.rate[, "y.f.1"],
        diag.subpop.rate[, "y.f.2"],
        diag.subpop.rate[, "y.f.3"],
        na.rm = TRUE
      )
    max.o.f <-
      100 * max(
        out.diag.o.f.1$value,
        out.diag.o.f.2$value,
        out.diag.o.f.3$value,
        diag.subpop.rate[, "o.f.1"],
        diag.subpop.rate[, "o.f.2"],
        diag.subpop.rate[, "o.f.3"],
        na.rm = TRUE
      )
    max.diag <- max(max.y.m, max.o.m, max.y.f, max.o.f) + 0.1
    max.msm <-
      max(out.diag.y.msm$value, out.diag.o.msm$value, na.rm = TRUE) * 100 + 0.1


    #get max values of incidence for setting y-axis in plots
    max.inc.y.m <-
      max(out.inc.y.m.1$value,
          out.inc.y.m.2$value,
          out.inc.y.m.3$value,
          na.rm = TRUE)
    max.inc.o.m <-
      max(out.inc.o.m.1$value,
          out.inc.o.m.2$value,
          out.inc.o.m.3$value,
          na.rm = TRUE)
    max.inc.y.f <-
      max(out.inc.y.f.1$value,
          out.inc.y.f.2$value,
          out.inc.y.f.3$value,
          na.rm = TRUE)
    max.inc.o.f <-
      max(out.inc.o.f.1$value,
          out.inc.o.f.2$value,
          out.inc.o.f.3$value,
          na.rm = TRUE)
    max.inc <-
      ceiling(max(max.inc.y.m, max.inc.o.m, max.inc.y.f, max.inc.o.f))
    max.inc.msm <-
      ceiling(max(out.inc.y.msm$value, out.inc.o.msm$value, na.rm = TRUE) + 1)



    out.inc.msm <- out.inc.y.msm
    out.inc.msm$value <-
      (out.inc.y.msm$value * sum(n.i[y.m.4])) + (out.inc.o.msm$value * sum(n.i[o.m.4]))
    out.inc.msm$value <- out.inc.msm$value / sum(n.i[m4])

    out.inc.msm %>% group_by(Var1) %>%
      summarize(mean = mean(value)) -> out.inc.msm.mean
  })

}

#' Layout the Target and Fit Plots
#'
#' Later I may turn this into a function that renders PDFs using these layouts.
#' Or I might return a list of page-by-page layouts for rendering in a PDF.
#' @import patchwork
#' @import lemon
layout_target_and_fit_plots <- function(l) {

  # (( l$plot.prev + grid_arrange_shared_legend(l$plot.age, l$plot.subpop) ) + plot_layout(ncol = 2, widths = c(1,2)) )+
  #   l$plot.p.symp + grid_arrange_shared_legend(l$plot.p.msm.y, l$plot.p.msm.o)

  plot_grid(
    l$plot.prev,
    grid_arrange_shared_legend(l$plot.age, l$plot.subpop),
    l$plot.p.symp,
    grid_arrange_shared_legend(l$plot.p.msm.y, l$plot.p.msm.o),
    ncol = 2,
    rel_widths = c(1.25, 2),
    labels = LETTERS[1:4]
  )

}
