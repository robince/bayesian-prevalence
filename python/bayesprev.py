import numpy as np
import scipy as sp
import scipy.stats
from scipy.stats import beta
from scipy.optimize import fsolve

# parameters for Gamma prior
r = 1
s = 1


def map(k, n, a=0.05, b=1):
    """ Bayesian maximum a posteriori estimate of population prevalence gamma
    under a uniform prior

    k : number of participants significant out of 
    n : total number of participants
    a : alpha value of within-participant test (default=0.05)
    b : sensitivity/beta of within-participant test (default=1)

    """

    theta = (k+r-1.0)/(n+r+s-2.0)
    if theta <= a:
        return 0.
    elif theta >= b:
        return 1.
    else:
        return (theta-a)/(b-a)

def posterior(x, k, n, a=0.05, b=1):
    """  Bayesian posterior of population prevalence gamma
    under a uniform prior

    x : values of gamma at which to evaluate the posterior density
    k : number of participants significant out of 
    n : total number of participants
    a : alpha value of within-participant test (default=0.05)
    b : sensitivity/beta of within-participant test (default=1)

    """

    theta = a+(b-a)*x
    post = (b-a)*beta.pdf(theta, k+r, n-k+s)
    post = post / (beta.cdf(b, k+r, n-k+s) - beta.cdf(a, k+r, n-k+s))
    return post

def bound(p, k, n, a=0.05, b=1):
    """Bayesian lower bound of population prevalence gamma under a uniform prior

    p : density the lower bound should bound (e.g. 0.95)
    k : number of participants significant out of 
    n : total number of participants
    a : alpha value of within-participant test (default=0.05)
    b : sensitivity/beta of within-participant test (default=1)

    """
    
    b1 = k+r
    b2 = n-k+s
    cdfp = (1-p)*beta.cdf(b,b1,b2) + p*beta.cdf(a,b1,b2)
    the_c = beta.ppf(cdfp,b1,b2)
    g_c = (the_c - a)/(b-a)
    return g_c

def hpdi(p, k, n, a=0.05, b=1):
    """ Bayesian highest posterior density interval of population prevalence gamma
    under a uniform prior

    p : HPDI to return (e.g. 0.95 for 95%)
    k : number of participants significant out of 
    n : total number of participants
    a : alpha value of within-participant test (default=0.05)
    b : sensitivity/beta of within-participant test (default=1)

    """

    b1 = k+r
    b2 = n-k+s

    # truncated beta pdf/cdf/icdf
    pdf = lambda x: beta.pdf(x,b1,b2) / (beta.cdf(b,b1,b2)-beta.cdf(a,b1,b2))
    cdf = lambda x: (beta.cdf(x,b1,b2)-beta.cdf(a,b1,b2)) / (beta.cdf(b,b1,b2)-beta.cdf(a,b1,b2))
    icdf = lambda x: beta.ppf( (1-x)*beta.cdf(a,b1,b2) + x*beta.cdf(b,b1,b2),b1,b2 )

    if k==a:
        x = np.array([a, icdf(p)])
    elif k==n:
        x = np.array([icdf(1-p), b])
    else:
        f = lambda x:  np.array([cdf(x[1])-cdf(x[0])-p, pdf(x[1])-pdf(x[0])]);
        x, info, ier, mesg = fsolve(f, np.array([icdf((1-p)/2), icdf((1+p)/2)]),full_output=True )

    # limit to valid theta values
    if (x[0]<a) or (x[1]<x[0]):
        x = np.array([a, icdf(p)])
    if x[1]>b: 
        x = np.array([icdf(1-p), b])
    hpdi = (x-a)/(b-a)
    return hpdi

