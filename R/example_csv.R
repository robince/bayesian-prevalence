  #
  # Example of loading data from a CSV and applying Bayesian prevalence


  # example_data.csv is a file with one binary value for each experimental
  # unit (participant, neuron, voxel etc.) where 1 indicates the within-unit
  # null hypothesis was rejected and 0 indicates it was not rejected
  
  # If the package "nleqslv" is not installed on your system, then
  # uncomment the following line and install this package.

  # install.packages("nleqslv", dependencies =TRUE)


   
 #    Load the nonlinear equation solver
   
 library(nleqslv)
         
 #    Set the working directory
 setwd("/users/jk/Desktop/PrevStuff")
         
 #  Load the required functions from the current working directory

 source("bayesprev.R")
   

  # load the data
  
  sigdat = read.csv("example_data.csv", header =FALSE);

  alpha = 0.05                       #   this specifies the alpha value used for the within-unit tests
  Ntests = nrow(sigdat);        #   number of tests (e.g. participants)
  Nsigtests = sum(sigdat);    #   number of significant tests
  
  # Plot the posterior  pdf for the population prevalence proportion
  
  xvals <- seq(0, 1, .01)
  pdf <- bayesprev_posterior(xvals, Nsigtests, Ntests)
  plot(xvals, pdf, type ="l", xlab = expression(gamma), ylab ="Posterior density", lwd=3)
  
  # 0.95 lower bound for the population prevalence proportion
  
  b = bayesprev_bound(0.05, Nsigtests, Ntests)
  print(b)
  
  # MAP estimate of the population prevalence proportion
  
  m = bayesprev_map(Nsigtests, Ntests)
  print(m)
  
  # 96% HPDI for the population prevalence proportion
  
  int = bayesprev_hpdi(0.96, Nsigtests, Ntests)
  print(int)
  
  
  
  # Example of possible Bayesian prevalence analyses on a common plot.
  

  xvals <- seq(0, 1, .01)
  pdf <- bayesprev_posterior(xvals, Nsigtests, Ntests)
  plot(xvals, pdf, type ="l", xlab = expression(gamma), ylab ="Posterior density", lwd=3)

  
  # Add the MAP estimate as a point
  
  xmap = bayesprev_map(Nsigtests, Ntests)
  pmap = bayesprev_posterior(xmap, Nsigtests, Ntests)
  points(xmap, pmap, cex =2, col= "red", pch=16)
  
  # Add the .95 lower bound as a vertical line
  
  xbound = bayesprev_bound(0.05, Nsigtests, Ntests)
  pbound = bayesprev_posterior(xbound, Nsigtests, Ntests)
  lines(c(xbound, xbound), c(0, pbound+0.5), col="blue", lwd=3)
  
  
  # Add the 0.96 HPDI
  
  int = bayesprev_hpdi(0.96, Nsigtests, Ntests)
  i1 = int[1]
  i2 = int[2]
  h1 = bayesprev_posterior(i1, Nsigtests, Ntests)
  h2 = bayesprev_posterior(i2, Nsigtests, Ntests)
  lines(c(i1, i2), c(h1, h2), col ="green", lwd=3)
  
  
  

 

  
  
  
  
