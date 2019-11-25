devtools::load_all()
library(here)
library(cowplot)

# load_start('SF')
# post.sample <- readRDS(here('inst/calibration_outcomes/SF_posterior_sample.rds'))
# pred_SF <- pred <- as.data.frame(post.sample$outputs)

render_4x4_fit_to_cases <- function(site) { 

load_start(site)
post.sample <- readRDS(here(paste0('inst/calibration_outcomes/', site, '_posterior_sample.rds')))
pred <- as.data.frame(post.sample$outputs)

attach(gc_env)

diag.age.sex.rate <- 
	if (!"year" %in% colnames(as.data.frame(diag.age.sex.rate))) {
		tibble::rownames_to_column(as.data.frame(diag.age.sex.rate), "year") 
	} else {
		diag.age.sex.rate
	}

diag.dat <- melt(diag.age.sex.rate[, 1:5], 
								 id.vars = "year", 
								 measured.vars = c("m_y", "m_o", "f_y", "f_o"))

diag.dat$cat <- 
	ifelse(diag.dat$variable == "m_y", "M 15-24 y", 
		 ifelse(diag.dat$variable == "m_o", "M 25-39 y", 
				ifelse(diag.dat$variable == "f_y", "F 15-24 y", "F 25-39 y")))


length.diag <- nrow(diag.age.sex.rate)
out.diag.y.m <-melt(as.matrix(subset(pred, select=prev.fit.diag.rate1:(prev.fit.diag.rate1+length.diag-1)))) # reported cases rates by age and sex only
out.diag.o.m <-melt(as.matrix(subset(pred, select=(prev.fit.diag.rate1+length.diag):(prev.fit.diag.rate1+2*length.diag-1))))
out.diag.y.f <-melt(as.matrix(subset(pred, select=(prev.fit.diag.rate1+2*length.diag):(prev.fit.diag.rate1+3*length.diag-1))))
out.diag.o.f<-melt(as.matrix(subset(pred, select=(prev.fit.diag.rate1+3*length.diag):(prev.fit.diag.rate1+4*length.diag-1))))

#get max values of reported cases for setting y-axis
max.y.m <- 100 * max(out.diag.y.m$value,
										 # diag.subpop.rate[,"y.m.1"],
										 # diag.subpop.rate[,"y.m.2"],
										 # diag.subpop.rate[,"y.m.3"],
										 na.rm = TRUE)

max.o.m <- 100 * max(out.diag.o.m$value,
										 # diag.subpop.rate[,"o.m.1"],
										 # diag.subpop.rate[,"o.m.2"],
										 # diag.subpop.rate[,"o.m.3"],
										 na.rm = TRUE)

max.y.f <- 100 * max(out.diag.y.f$value,
										 # diag.subpop.rate[,"y.f.1"],
										 # diag.subpop.rate[,"y.f.2"],
										 # diag.subpop.rate[,"y.f.3"],
										 na.rm = TRUE)

max.o.f <- 100 * max(out.diag.o.f$value,
										 # diag.subpop.rate[,"o.f.1"],
										 # diag.subpop.rate[,"o.f.2"],
										 # diag.subpop.rate[,"o.f.3"],
										 na.rm = TRUE)

# max.diag <- max(max.y.m, max.o.m, max.y.f, max.o.f) + .1
max.diag <- 2.7

diag.data <- cbind(
	year=(seq(as.numeric(diag.age.sex.rate[1,"year"]),as.numeric(diag.age.sex.rate[nrow(diag.age.sex.rate), "year"]),1)), 
	as.data.frame(diag.rate), 
	cat=rep(c("y.m", "o.m", "y.f", "o.f"), each=nrow(diag.age.sex.rate)))


# Plot Reported & Diagnosed Cases and Corresponding Model Fits for Young Men
plot.diag.y.m <- ggplot(data=out.diag.y.m)+
	geom_line(aes(x=Var2, y=value*100, group=Var1), color="grey") +
	stat_summary(aes(x=Var2, y=value*100, group=1), geom="line", fun.y="median", color="black", size=1) +
	theme_classic() +
	theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
	geom_point(data=diag.data[1:length.diag,],aes(x=1:length.diag, y=diag.rate*100), color="red", shape=15, size=2) +
	labs(title="All M: 15-24y", x="Year", y="Reported diagnoses per 100") +
	coord_cartesian(ylim=c(0,max.diag))+
	scale_x_discrete(labels=seq((end.year-length.diag+1),end.year,1))

# Plot Reported & Diagnosed Cases and Corresponding Model Fits for Old Men
plot.diag.o.m <- ggplot(data=out.diag.o.m)+
	geom_line(aes(x=Var2, y=value*100, group=Var1), color="grey") +
	stat_summary(aes(x=Var2, y=value*100, group=1), geom="line", fun.y="median", color="black", size=1) +
	theme_classic() +
	theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
	geom_point(data=diag.data[(length.diag+1):(2*length.diag),],aes(x=1:length.diag, y=diag.rate*100), color="red", shape=15, size=2) +
	labs(title="All M: 25-39y", x="Year", y="Reported diagnoses per 100") +
	coord_cartesian(ylim=c(0,max.diag))+
	scale_x_discrete(labels=seq((end.year-length.diag+1),end.year,1))

# Plot Reported & Diagnosed Cases and Corresponding Model Fits for Young Women 
plot.diag.y.f <- ggplot(data=out.diag.y.f)+
	geom_line(aes(x=Var2, y=value*100, group=Var1), color="grey") +
	stat_summary(aes(x=Var2, y=value*100, group=1), geom="line", fun.y="median", color="black", size=1) +
	theme_classic() +
	theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
	geom_point(data=diag.data[(2*length.diag+1):(3*length.diag),],aes(x=1:length.diag, y=diag.rate*100), color="red", shape=15, size=2) +
	labs(title="All F: 15-24y", x="Year", y="Reported diagnoses per 100") +
	coord_cartesian(ylim=c(0,max.diag))+
	scale_x_discrete(labels=seq((end.year-length.diag+1),end.year,1))

# Plot Reported & Diagnosed Cases and Corresponding Model Fits for Old Women 
plot.diag.o.f <- ggplot(data=out.diag.o.f)+
	geom_line(aes(x=Var2, y=value*100, group=Var1), color="grey") +
	stat_summary(aes(x=Var2, y=value*100, group=1), geom="line", fun.y="median", color="black", size=1) +
	theme_classic() +
	theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) +
	geom_point(data=diag.data[(3*length.diag+1):(4*length.diag),],aes(x=1:length.diag, y=diag.rate*100), color="red", shape=15, size=2) +
	labs(title="All F: 25-39y", x="Year", y="Reported diagnoses per 100") +
	coord_cartesian(ylim=c(0,max.diag))+
	scale_x_discrete(labels=seq((end.year-length.diag+1),end.year,1))


p <- plot_grid(plot.diag.y.m, plot.diag.o.m, plot.diag.y.f, plot.diag.o.f)
sitenames <- c(BA = "Baltimore", SF = 'San Francisco')
ordering_of_sites <- c(BA = 'A)', SF = 'B)')
title <- ggdraw() + draw_label(paste0(ordering_of_sites[[site]], " Reported Diagnoses and Model Fit in ", sitenames[[site]]), fontface = 'bold', size = 12)
p1 <- plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins

ggsave(p1, filename = paste0("4x4_fit_to_cases_", site, ".png"), width = unit(4.5, 'in'), height = unit(4, 'in'))
detach(gc_env)
}

render_4x4_fit_to_cases('BA')
render_4x4_fit_to_cases('SF')

