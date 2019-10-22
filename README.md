# bayesian-prevalence

Bayesian inference of population prevalence.

This package includes code for Matlab, Python and R implementing Bayesian prevalence inference. 

Consider that a statistical test (with common false positive rate alpha) is performed in each participant in a psychology or neuroimaging experiment, or on each single unit recorded in an electrophysiology experiment. Following this first-level analysis, we can compute the Bayesian estimate of the prevalence of the effect in the population using only three numbers: the total number of tests `n`, out of which `k` are positive, with false positive rate (alpha) `a`. These numbers can be specified directly in a function call or obtained from a variable that indicates the result of applying the individual tests at the first level. The `example_csv` scripts give an example of loading this first-level within-participant signficance data and applying the second-level prevalence functions. 
 
`bayesprev_example.{m,R,py}` simulates data under a hierarchical normal model, applies a t-test against zero within each participant at the first level, and applies Bayesian prevalence inference at the second level. Users can adapt this example to load their own raw data, replace the t-test with any other within-participant statistical test, or load the indicator variable for significance directly and apply the second level test (see also `example_csv`).
 
These functions assume a uniform prior on the prevalence proportion. The inferred prevalence value gives an estimate of the proportion of the population that has the tested effect. For mathematical details of the Bayesian inference procedure see BayesianPrevalence.pdf 

## Functions

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
Gives the posterior density of the population prevalence proportion computed at values in `x` (vector with values `0<x<1`). E.g. (Matlab) `x=linspace(0,1,100);plot(x,bayesprev_posterior(x,k,n,a))`

