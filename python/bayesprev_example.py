# Example of how to use Bayesian prevalence functions
#
# 1.   Simulate or load within-participant raw experimental data
# 2.   LEVEL 1: Apply statstical test at the individual level
# 3.   LEVEL 2: Apply Bayesian Prevalence to the outcomes of Level 1

import numpy as np
import scipy as sp
import bayesprev
import matplotlib.pyplot as plt

#
# 1.   Simulate or load within-participant raw experimental data
#

# 1.1. Simulate within-participant raw experimental data
Nsub = 20 # number of particpants
Nsamp = 100 # trials/samples per participant
sigma_w = 10 # within-participant SD
sigma_b = 2 # between-participant SD
mu_g = 1 # population mean

# per participant mean drawn from population normal distribution
submeanstrue = np.random.normal(mu_g, sigma_b, Nsub)
# rawdat holds trial data for each participant 
rawdat = np.zeros((Nsamp, Nsub))
for si in range(Nsub):
    # generate trials for each participant
    rawdat[:,si] = np.random.normal(submeanstrue[si], sigma_w, Nsamp)

# 1.2.Load within-participant raw experimental data
# Load your own data into the variable rawdat with dimensions [Nsamp Nsub],
# setting Nsamp and Nsub accordingly.

#
# 2.  LEVEL 1
#

# 2.1. Within-participant statistical test
# This function performs within-participant statistical test. Here, a t-test for
# non-zero mean which is the simplest statistical test. In general, any
# statistical test can be used at Level 1.

# calculates a t-test against 0 mean independently for each participant
[t, p] = sp.stats.ttest_1samp(rawdat,0)
# p holds p-values of test for each participant
alpha = 0.05 # false positive rate of test
indsig = p<alpha
# the binary variable indsig indicates whether the within-participant
# t-test is, or not, significant for each participant (1 entry for each
# participant)

# 2.2. Loading within-participant statistical test.
# You can also load your own within-participant statistical test results here.
# Load binary results into binary indsig vector with one entry per participant.
# See also example_csv.py for an example.

#
# 3.  LEVEL 2
#

# Bayesian prevalence inference is performed with three numbers: 
# k, the number of significant participants (e.g. sum of binary indicator
# variable)
# n, the number of participants in the sample
# alpha, the false positive rate

n = Nsub
k = indsig.sum()
# alpha defined above


f = plt.figure()
col = '#1f77b4'

# plot posterior distribution of population prevalence
x = np.linspace(0,1,100)
posterior = bayesprev.posterior(x,k,n,alpha)
plt.plot(x,posterior,color=col)

# add MAP as a point
xmap = bayesprev.map(k,n,alpha)
pmap = bayesprev.posterior(xmap,k,n,alpha)
plt.plot(xmap, pmap,'.',color=col,markersize=20)

# add lower bound as a vertical line
bound = bayesprev.bound(0.95,k,n,alpha)
pbound = bayesprev.posterior(bound,k,n,alpha)
plt.vlines(bound, 0, pbound, linestyles=':',colors=col) 

# add 95% HPDI
oil = 2;
iil = 4;
h = bayesprev.hpdi(0.95,k,n,alpha)
plt.hlines(pmap, h[0], h[1],linewidths=oil,colors=col)
# add 50% HPDI
h = bayesprev.hpdi(0.5,k,n,alpha)
plt.hlines(pmap, h[0], h[1],linewidths=iil,colors=col)

plt.xlabel('Population prevalence proportion')
plt.ylabel('Posterior density')

plt.show()
