# Example of how to use Bayesian prevalence functions
#
# 1.   Simulate or load within-participant raw experimental data
# 2.   LEVEL 1: Apply statstical test at the individual level
# 3.   LEVEL 2: Apply Bayesian Prevalence to the outcomes of Level 1

source("bayesprev.R")

#
# 1.   Simulate or load within-participant raw experimental data
#

# 1.1. Simulate within-participant raw experimental data
Nsub<-20 # number of particpants
Nsamp <- 100 # trials/samples per participant
sigma_w <- 10 # within-participant SD
sigma_b <- 2 # between-participant SD
mu_g <- 1 # population mean

# per participant mean drawn from population normal distribution
submeanstrue <- rnorm(Nsub, mu_g, sigma_b)
# rawdat holds trial data for each participant 
rawdat = matrix(0, nrow=Nsamp, ncol=Nsub)
for( si in 1:Nsub ){
  # generate trials for each participant
  rawdat[, si] = rnorm(Nsamp, submeanstrue[si], sigma_w)
}

# 1.2.Load within-participant raw experimental data
# Load your own data into the variable rawdat with dimensions [Nsamp Nsub],
# setting Nsamp and Nsub accordingly.

#
# 2.  LEVEL 1
#

# 2.1. Within-participant statistical test
# This loop performs within-participant statistical test. Here, a t-test for
# non-zero mean which is the simplest statistical test. In general, any
# statistical test can be used at Level 1.

p <- vector(mode="numeric",length=Nsub)
for( si in 1:Nsub ){
  t = t.test(rawdat[,si], mu=0)
  p[si] = t[3]
}
# p holds p-values of test for each participant
alpha = 0.05  
indsig = p<alpha
# the binary variable indsig indicates whether the within-participant
# t-test is, or not, significant for each participant (1 entry for each
# participant)

# 2.2. Loading within-participant statistical test.
# You can also load your own within-participant statistical test results here.
# Load binary results into binary indsig vector with one entry per participant. 
# See also example_csv.R for an example.

#
# 3.  LEVEL 2
#

# Bayesian prevalence inference is performed with three numbers: 
# k, the number of significant participants (e.g. sum of binary indicator
# variable)
# n, the number of participants in the sample
# alpha, the false positive rate
n <- Nsub
k <- sum(indsig)

# plot posterior distribution of population prevalence
xvals <- seq(0, 1, .01)
pdf <- bayesprev_posterior(xvals, k, n)
plot(xvals, pdf, type ="l", xlab = expression(gamma), ylab ="Posterior density", lwd=3)

# Add the MAP estimate as a point
xmap = bayesprev_map(k, n)
pmap = bayesprev_posterior(xmap, k, n)
points(xmap, pmap, cex =2, col= "red", pch=16)

# Add the .95 lower bound as a vertical line
xbound = bayesprev_bound(0.95, k, n)
pbound = bayesprev_posterior(xbound, k, n)
lines(c(xbound, xbound), c(0, pbound+0.5), col="blue", lwd=3)

# Add the 0.96 HPDI
int = bayesprev_hpdi(0.96, k, n)
i1 = int[1]
i2 = int[2]
h1 = bayesprev_posterior(i1, k, n)
h2 = bayesprev_posterior(i2, k, n)
lines(c(i1, i2), c(h1, h2), col ="green", lwd=3)



