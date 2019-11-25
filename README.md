
# Modeling Gonorrhea Transmission in Baltimore and San Francisco

The `gcRegional` R package implemented here defines a model and analyzes outcomes
from this model for the purpose of providing simulation model evidence for
the forthcoming research manuscript *The potential population-level impact of
different gonorrhea screening strategies in Baltimore and San Francisco:  an
exploratory mathematical modeling analysis*. 

# Authorship 

The author list includes: 

> Minttu M. Rönn<sup>#*1</sup>, Christian Testa<sup>#1</sup>, Ashleigh R.
> Tuite<sup>1</sup>, Harrell W. Chesson<sup>2</sup>, Thomas L. Gift<sup>2</sup>,
> Christina Schumacher<sup>3,4</sup>, Sarah L. Williford<sup>3,4</sup>, Lin
> Zhu<sup>1</sup>, Meghan Bellerose<sup>1</sup>, Rebecca Earnest<sup>1</sup>,
> Yelena Malyuta<sup>1</sup>, Katherine K. Hsu<sup>5</sup>, Joshua
> A. Salomon<sup>6</sup>, Nicolas A. Menzies<sup>1

<sup>\#</sup> Contributed equally <br>
<sup>\*</sup> Corresponding author

Affiliations:

1.	Department of Global Health, Harvard T. H. Chan School of Public Health
2.	Centers for Disease Control and Prevention, Division of STD Prevention
3.	Johns Hopkins University School of Medicine
4.	Baltimore City Health Department 
5.	Division of STD Prevention & HIV/AIDS Surveillance, Massachusetts Department of Public Health
6.	Center for Health Policy / Center for Primary Care and Outcomes Research, Stanford University

# Funding Source

This work was funded by the U.S. Centers for Disease Control and Prevention,
National Center for HIV, Viral Hepatitis, STD, and TB Prevention Epidemiologic
and Economic Modeling Agreement (NEEMA, #5U38PS004644)

# Data Sources

In this project we used data from the NHANES, NSFG, and SSuN data surveys.

# Methodology

The `gcRegional` R package implements a deterministic compartmental model 
of gonorrhea transmission. 

The approach taken here is to write the simulation module (i.e. a differential
equations solver) in C++ accessible through an Rcpp interface in this package.
This simulation model is calibrated to the observed data, and using calibrated
models we are able to estimated the impact of changes in screening patterns on
the fit models' future outcomes. 

We take a Bayesian modeling approach to calibrate our model; We formally construct
the dLogPosterior function and estimate the parametric uncertainty of our 
model using an adaptive MCMC algorithm. 

# Installing the Package

### Installation Prerequisites

One must have Rcpp installed and configured correctly in order to install and use `gcRegional`. 

### Install gcRegional 

Installing this package will allow the user to run the model and to produce 
plots such as those we have included in the manuscript and in this package. 

Install the package and its dependencies using `library(devtools);
install_github("PPML/gc-regional.git", dependencies=T)`.

### Cloning the Project for Calibration

Fitting the model and storing the outcomes of the calibration process is
functionality that we advice one does by cloning this git repository into
a local directory. 

This is because the code to store calibrated results uses the `here` package
in the context of the `gcRegional` package to store results in the `inst/`
directory to make them available in the package. 

# Running the Model

```
library(gcRegional) # load our package

# Use the load_start function to load the model configurations for either
# Baltimore or San Francisco.
load_start('BA') # or 'SF'

# The parametrization for the model in the selected location is given by 
# parameters in the environment `gc_env`.

# > str(gc_env$theta)
#  Named num [1:106] -0.0306 -0.7014 -0.0688 -0.9072 -0.1881 ...
#  - attr(*, "names")= chr [1:106] "logit.epsilon.1" "logit.epsilon.2" "logit.epsilon.3" "logit.epsilon.4" ...

e <- create_gcSim_param_env(theta=gc_env$theta)
run_gcSim_with_environment(e)
```

