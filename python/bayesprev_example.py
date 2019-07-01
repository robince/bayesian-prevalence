# Example of how to use Bayesian prevalence functions

import numpy as np
import scipy as sp
import bayesprev
import matplotlib.pyplot as plt


# simulate some data
Nsub = 20 # number of particpants
Nsamp = 100 # trials/samples per participant
sigma_w = 10 # within-participant SD
sigma_b = 2 # between-participant SD
mu_g = 1 # population mean

# generate within-participant data
submeanstrue = np.random.normal(mu_g, sigma_b, Nsub)
rawdat = np.zeros((Nsamp, Nsub))
for si in range(Nsub):
    rawdat[:,si] = np.random.normal(submeanstrue[si], sigma_w, Nsamp)

# perform within-participant statistical test (here t-test)
[t, p] = sp.stats.ttest_1samp(rawdat,0)

# within-participant test results for prevalence inference 
alpha = 0.05  
n = Nsub
indsig = p<alpha
k = indsig.sum()

# plot posterior distribution of population prevalence

f = plt.figure()
col = '#1f77b4'

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
h = bayesprev.hpdi(0.5,k,n,alpha)
plt.hlines(pmap, h[0], h[1],linewidths=iil,colors=col)

plt.xlabel('Population prevalence proportion')
plt.ylabel('Posterior density')

plt.show()
