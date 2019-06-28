# bayesian-prevalence

Bayesian inference of population prevalence.

This package includes code for Matlab, Python and R implementing Bayesian prevalence inference. 

This can be applied whenever a number of statistical tests are performed with a common false positive rate alpha. For example, a test in each participant during a psychology or neuroimaging experiment, or a test on each single unit recorded during an electrophysiology experiment.

The total number of tests is `n`, out of which `k` are positive and each test has false positive rate (alpha) `a`. 

These functions assume a uniform prior on the prevalence proportion. The inferred prevalence value gives an estimate of the proportion of the population which have the tested effect (i.e. would be true positives if tested).  

Then:

`bayesprev_map(k,n,a)` *Matlab*, *R*  
`bayesprev.map(k,n,a)` *Python*  
Gives the maximum a posteriori estimate of the population prevalence proportion (the mode or most likely value of the posterior distribution).

`bayesprev_bound(p,k,n,a)` *Matlab*, *R*  
`bayesprev.bound(p,k,n,a)` *Python*  
Gives a lower bound p for the population prevalence proportion from the `1-p` quantile of the posterior distribution. (`0<p<1`)

`bayesprev_hpdi(p,k,n,a)` *Matlab*, *R*  
`bayesprev.hpdi(p,k,n,a)` *Python*  
Gives the highest posterior density interval for the population prevalence proportion. This is the shortest interval of prevalence proportions for which the integral of the posterior density is `p`. (`0<p<1`).

`bayesprev_posterior(x,k,n,a)` *Matlab*, *R*  
`bayesprev.posterior(x,k,n,a)` *Python*  
Gives the posterior density of the population prevalence proportion computed at values in `x` (vector with values `0<x<1`).


