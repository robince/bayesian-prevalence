"""Bayesian estimation of population prevalence

Bayesian inference of population prevalence
RAA Ince, AT Paton, JW Kay & PG Schyns
(2021) eLife 10:e62461 doi: 10.7554/eLife.62461

If a statistical test with false positive rate alpha is performed in 
N participants (or other replication units) and k are found to be significant
then we can estimate the population prevalence of a true positive effect. 
This is the population level within-participant replication probability. 

"""
__version__ = "0.1.1"

import numpy as np
import scipy as sp
from scipy.optimize import fsolve
from scipy.stats import beta

# parameters for Gamma prior
r = 1
s = 1


def map(k, n, a=0.05, b=1):
    """Bayesian maximum a posteriori estimate of population prevalence gamma
    under a uniform prior

    k : number of participants significant out of
    n : total number of participants
    a : alpha value of within-participant test (default=0.05)
    b : sensitivity/beta of within-participant test (default=1)

    """

    theta = (k + r - 1.0) / (n + r + s - 2.0)
    if theta <= a:
        return 0.0
    elif theta >= b:
        return 1.0
    else:
        return (theta - a) / (b - a)


def posterior(x, k, n, a=0.05, b=1):
    """Bayesian posterior of population prevalence gamma
    under a uniform prior

    x : values of gamma at which to evaluate the posterior density
    k : number of participants significant out of
    n : total number of participants
    a : alpha value of within-participant test (default=0.05)
    b : sensitivity/beta of within-participant test (default=1)

    """

    theta = a + (b - a) * x
    post = (b - a) * beta.pdf(theta, k + r, n - k + s)
    post = post / (beta.cdf(b, k + r, n - k + s) - beta.cdf(a, k + r, n - k + s))
    return post


def bound(p, k, n, a=0.05, b=1):
    """Bayesian lower bound of population prevalence gamma under a uniform prior

    p : density the lower bound should bound (e.g. 0.95)
    k : number of participants significant out of
    n : total number of participants
    a : alpha value of within-participant test (default=0.05)
    b : sensitivity/beta of within-participant test (default=1)

    """

    b1 = k + r
    b2 = n - k + s
    cdfp = (1 - p) * beta.cdf(b, b1, b2) + p * beta.cdf(a, b1, b2)
    the_c = beta.ppf(cdfp, b1, b2)
    g_c = (the_c - a) / (b - a)
    return g_c


def hpdi(p, k, n, a=0.05, b=1):
    """Bayesian highest posterior density interval of population prevalence gamma
    under a uniform prior

    p : HPDI to return (e.g. 0.95 for 95%)
    k : number of participants significant out of
    n : total number of participants
    a : alpha value of within-participant test (default=0.05)
    b : sensitivity/beta of within-participant test (default=1)

    """

    b1 = k + r
    b2 = n - k + s

    # truncated beta pdf/cdf/icdf
    tbpdf = lambda x: beta.pdf(x, b1, b2) / (beta.cdf(b, b1, b2) - beta.cdf(a, b1, b2))
    tbcdf = lambda x: (beta.cdf(x, b1, b2) - beta.cdf(a, b1, b2)) / (
        beta.cdf(b, b1, b2) - beta.cdf(a, b1, b2)
    )
    tbicdf = lambda x: beta.ppf((1 - x) * beta.cdf(a, b1, b2) + x * beta.cdf(b, b1, b2), b1, b2)

    if k == a:
        x = np.array([a, tbicdf(p)])
    elif k == n:
        x = np.array([tbicdf(1 - p), b])
    else:
        f = lambda x: np.array([tbcdf(x[1]) - tbcdf(x[0]) - p, tbpdf(x[1]) - tbpdf(x[0])])
        x, info, ier, mesg = fsolve(
            f, np.array([tbicdf((1 - p) / 2), tbicdf((1 + p) / 2)]), full_output=True
        )

    # limit to valid theta values
    if (x[0] < a) or (x[1] < x[0]):
        x = np.array([a, tbicdf(p)])
    if x[1] > b:
        x = np.array([tbicdf(1 - p), b])
    hpdi = (x - a) / (b - a)
    return hpdi


def diff_between(k1, n1, k2, n2, p=0.96, a=0.05, b=1, Nsamp=10000):
    """Bayesian maximum a posteriori estimate of the difference in prevalence
    when the same test is applied to two groups

    k1 : number of participants significant in group 1 out of
    n1 : total number of participants in group 1
    k2 : number of participants significant in group 2 out of
    n2 : total number of participants in group 2
    p  : coverage for highest-posterior density interval (in [0 1])
    a  : alpha value of within-participant test (default=0.05)
    b  : sensitivity/beta of within-participant test (default=1)
    Nsamp : number of samples from the posterior

    Outputs:
    map    : maximum a posteriori estimate of the difference in prevalence:
             gamma_1 - gamma_2
    post_x : x-axis for kernel density fit of posterior distribution of the
             above
    post   : posterior distribution from kernel density fit
    hpdi   : highest-posterior density interval with coverage p
    probGT : estimated posterior probability that the prevalence is higher in group 1
    logoddsGT : estimated log odds in favour of the hypothesis that the prevalence is higher in group 1
    samples : posterior samples

    """

    # gamma priors = Beta(r,s)
    r1 = 1
    s1 = 1
    r2 = 1
    s2 = 1

    # Parameters for Beta posteriors
    m11 = k1 + r1
    m12 = n1 - k1 + s1
    m21 = k2 + r2
    m22 = n2 - k2 + s2

    # Generate truncated beta samples
    # fix numerical issue
    r1 = (beta.cdf(a, m11, m12), beta.cdf(b, m11, m12))
    r2 = (beta.cdf(a, m21, m22), beta.cdf(b, m21, m22))
    if np.any([np.isclose(*r, rtol=1e-12, atol=1e-12) for r in [r1, r2]]):
        res = {
            x: np.NaN for x in ["map", "post_x", "post", "hpdi", "probGT", "logoddsGT", "samples"]
        }
        return res

    th1 = beta.ppf(np.random.uniform(r1[0], r1[1], Nsamp), m11, m12)
    th2 = beta.ppf(np.random.uniform(r2[0], r2[1], Nsamp), m21, m22)

    # vector of estimates of prevalence differences
    samples = (th1 - th2) / (b - a)

    # kernel density estimate of posterior
    post_x = np.linspace(-1, 1, 200)
    kde = sp.stats.gaussian_kde(samples)
    post = kde(post_x)
    map = post_x[np.argmax(post)]

    # Estimate the posterior probability, and logodds, that the prevalence is higher for group 1.
    # Laplace's rule of succession used to avoid estimates of 0 or 1
    probGT = (np.sum(samples > 0) + 1) / (Nsamp + 2)
    logoddsGT = np.log(probGT / (1 - probGT))
    hpdi = _hpdi(samples, p)

    res = {
        "map": map,
        "post_x": post_x,
        "post": post,
        "hpdi": hpdi,
        "probGT": probGT,
        "logoddsGT": logoddsGT,
        "samples": samples,
    }
    return res


def diff_within(k11, k10, k01, n, p=0.96, a=0.05, b=1, Nsamp=10000):
    """Bayesian maximum a posteriori estimate of the difference in prevalence
    when the same test is applied to two groups

    k1 : number of participants significant in group 1 out of
    n1 : total number of participants in group 1
    k2 : number of participants significant in group 2 out of
    n2 : total number of participants in group 2
    p  : coverage for highest-posterior density interval (in [0 1])
    a  : alpha value of within-participant test (default=0.05)
    b  : sensitivity/beta of within-participant test (default=1)
    Nsamp : number of samples from the posterior

    Outputs:
    map    : maximum a posteriori estimate of the difference in prevalence:
             gamma_1 - gamma_2
    post_x : x-axis for kernel density fit of posterior distribution of the
             above
    post   : posterior distribution from kernel density fit
    hpdi   : highest-posterior density interval with coverage p
    probGT : estimated posterior probability that the prevalence is higher in group 1
    logoddsGT : estimated log odds in favour of the hypothesis that the prevalence is higher in group 1
    samples : posterior samples

    """

    # Parameters for the Dirichlet prior distribution (1,1,1,1) = uniform
    r11 = 1
    r10 = 1
    r01 = 1
    r00 = 1

    # parameters for posterior Dirichlet distribution
    k00 = n - k11 - k10 - k01
    m11 = k11 + r11
    m10 = k10 + r10
    m01 = k01 + r01
    m00 = k00 + r00

    r11 = (beta.cdf(0, m11, m10 + m01 + m00), beta.cdf(b, m11, m10 + m01 + m00))
    if np.any(np.isclose(*r11, rtol=1e-12, atol=1e-12)):
        res = {
            x: np.NaN for x in ["map", "post_x", "post", "hpdi", "probGT", "logoddsGT", "samples"]
        }
        return res
    # samples from the truncated Dirichlet posterior
    z11 = np.random.uniform(r11[0], r11[1], Nsamp)
    th11 = beta.ppf(z11, m11, m10 + m01 + m00)

    lo = np.maximum((a - th11) / (1 - th11), 0)
    hi = (b - th11) / (1 - th11)

    r10 = beta.cdf(lo, m10, m01 + m00), beta.cdf(hi, m10, m01 + m00)
    if np.any(np.isclose(*r10, rtol=1e-12, atol=1e-12)):
        res = {
            x: np.NaN for x in ["map", "post_x", "post", "hpdi", "probGT", "logoddsGT", "samples"]
        }
        return res
    z10 = np.random.uniform(r10[0], r10[1], Nsamp)
    u10 = beta.ppf(z10, m10, m01 + m00)
    th10 = (1 - th11) * u10

    lo = np.maximum((a - th11) / (1 - th11 - th10), 0)
    hi = np.minimum((b - th11) / (1 - th11 - th10), 1)
    r01 = (beta.cdf(lo, m01, m00), beta.cdf(hi, m01, m00))
    if np.any(np.isclose(*r01, rtol=1e-12, atol=1e-12)):
        res = {
            x: np.NaN for x in ["map", "post_x", "post", "hpdi", "probGT", "logoddsGT", "samples"]
        }
        return res
    z01 = np.random.uniform(r01[0], r01[1], Nsamp)
    u01 = beta.ppf(z01, m01, m00)
    th01 = (1 - th11 - th10) * u01

    th00 = 1 - th11 - th10 - th01

    # samples of posterior prevalence difference
    samples = (th10 - th01) / (b - a)

    # kernel density estimate of posterior
    post_x = np.linspace(-1, 1, 200)
    kde = sp.stats.gaussian_kde(samples)
    post = kde(post_x)
    map = post_x[np.argmax(post)]

    # Estimate the posterior probability, and logodds, that the prevalence is higher for group 1.
    # Laplace's rule of succession used to avoid estimates of 0 or 1
    probGT = (np.sum(samples > 0) + 1) / (Nsamp + 2)
    logoddsGT = np.log(probGT / (1 - probGT))
    hpdi = _hpdi(samples, p)

    res = {
        "map": map,
        "post_x": post_x,
        "post": post,
        "hpdi": hpdi,
        "probGT": probGT,
        "logoddsGT": logoddsGT,
        "samples": samples,
    }
    return res


def _hpdi(data, p):
    """HPDI modified from https://arviz-devs.github.io/arviz/"""
    data = data.flatten()
    n = len(data)
    data = np.sort(data)
    interval_idx_inc = int(np.floor(p * n))
    n_intervals = n - interval_idx_inc
    interval_width = data[interval_idx_inc:] - data[:n_intervals]

    if len(interval_width) == 0:
        raise ValueError("Too few elements for interval calculation.")

    min_idx = np.argmin(interval_width)
    hdi_min = data[min_idx]
    hdi_max = data[min_idx + interval_idx_inc]

    hpdi = np.array([hdi_min, hdi_max])

    return hpdi
