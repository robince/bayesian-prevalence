# 
# Example of loading data from a CSV and applying Bayesian prevalence
# second level

# example_data.csv is a file with one binary value for each experimental
# unit (participant, neuron, voxel etc.) where 1 indicates the within-unit
# null hypothesis was rejected and 0 indicates it was not rejected

import numpy as np
import scipy as sp
import bayesprev

# load the data
sigdat = np.loadtxt('example_data.csv')

# input values for Bayesian prevalence
alpha = 0.05 # this specifies the alpha value used for the within-unit tests
Ntests = sigdat.size # number of tests (e.g. participants)
Nsigtests = sigdat.sum() # number of significant tests

# example of Bayesian prevalence inference on a single plot
f = plt.figure()
k = Nsigtests
n = Ntests
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
# add 50% HPDI
h = bayesprev.hpdi(0.5,k,n,alpha)
plt.hlines(pmap, h[0], h[1],linewidths=iil,colors=col)

plt.xlabel('Population prevalence proportion')
plt.ylabel('Posterior density')

plt.show()
