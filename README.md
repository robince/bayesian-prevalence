# bayesian-prevalence

Bayesian inference of population prevalence.

This package includes code for Matlab, Python and R implementing Bayesian prevalence inference as described in:

Bayesian inference of population prevalence  
RAA Ince, AT Paton, JW Kay & PG Schyns  
(2021) *eLife* 10:e62461 [doi: 10.7554/eLife.62461](https://doi.org/10.7554/eLife.62461)

Consider that a statistical test (with common false positive rate alpha) is performed in each participant in a psychology or neuroimaging experiment, or on each single unit recorded in an electrophysiology experiment. Following this first-level analysis, we can compute the Bayesian estimate of the prevalence of true positive results to such a test in the population using only three numbers: the total number of tests `n`, out of which `k` are positive, with false positive rate (alpha) `a`. These numbers can be specified directly in a function call or obtained from a variable that indicates the result of applying the individual tests at the first level. The `example_csv` scripts give an example of loading this first-level within-participant signficance data and applying the second-level prevalence functions. 
 
`bayesprev_example.{m,R,py}` simulates data under a hierarchical normal model, applies a t-test against zero within each participant at the first level, and applies Bayesian prevalence inference at the second level. Users can adapt this example to load their own raw data, replace the t-test with any other within-participant statistical test, or load the indicator variable for significance directly and apply the second level test (see also `example_csv`).
 
These functions assume a uniform prior on the prevalence proportion. The inferred prevalence value gives an estimate of the proportion of the population that has the tested effect. For mathematical details of the Bayesian inference procedure see BayesianPrevalence.pdf 

There are also functions which return posterior estimates and samples for the difference in prevalence between two populations, or between two tests applied to the same sample. 

BayesianPrevalence.pdf and BayesianPrevalenceDifference.pdf give provide more mathematical detail of the methods. 

The folder `paper` contains Matlab scripts to generate the figures in the manuscript. 

## Functions

`bayesprev_map(k,n,a,b)` *Matlab*, *R*  
`bayesprev.map(k,n,a,b)` *Python*  
Gives the maximum a posteriori estimate of the population prevalence proportion (the mode or most likely value of the posterior distribution).

`bayesprev_bound(p,k,n,a,b)` *Matlab*, *R*  
`bayesprev.bound(p,k,n,a,b)` *Python*  
Gives a lower bound p for the population prevalence proportion from the `1-p` quantile of the posterior distribution. (`0<p<1`)

`bayesprev_hpdi(p,k,n,a,b)` *Matlab*, *R*  
`bayesprev.hpdi(p,k,n,a,b)` *Python*  
Gives the highest posterior density interval for the population prevalence proportion. This is the shortest interval of prevalence proportions for which the integral of the posterior density is `p`. (`0<p<1`).

`bayesprev_posterior(x,k,n,a,b)` *Matlab*, *R*  
`bayesprev.posterior(x,k,n,a,b)` *Python*  
Gives the posterior density of the population prevalence proportion computed at values in `x` (vector with values `0<x<1`). E.g. (Matlab) `x=linspace(0,1,100);plot(x,bayesprev_posterior(x,k,n,a))`

`bayesprev_diff_between(k1,n1,k2,n2,p,a,b,Nsamp)` *Matlab*, *R*  
`bayesprev.diff_between(k1,n1,k2,n2,p,a,b,Nsamp)` *Python*  
Returns MAP, kernel density fit to posterior density, HPDI of p, probability `g1-g2>0`, log odds `g1>g2`, and posterior samples for the population prevalence difference `g1-g2` between two populations. 

`bayesprev_diff_within(k11,k10,k01,n,p,a,b,Nsamp)` *Matlab*, *R*  
`bayesprev.diff_within(k11,k10,k01,n,p,a,b,Nsamp)` *Python*  
Returns MAP, kernel density fit to posterior density, HPDI of p, probability `g1-g2>0`, log odds `g1>g2`, and posterior samples for the population prevalence difference `g1-g2` for two tests applied to the same sample. 

## Installation

### Python

`pip install bayesprev`

### Matlab

Clone or download this repository, then in Matlab add folder to path: 

`addpath('/path/to/bayesprev/matlab')`

### R

Clone or download this repository, then in R:

`source("bayesprev.R")`

## Contact

Robin Ince, robince@gmail.com (Python code, MATLAB code)
Jim Kay, jimkay049@gmail.com (R code, Technical note on Bayesian prevalence)
