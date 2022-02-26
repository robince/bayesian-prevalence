
# nleqslv package required for HPDI optimization
library(nleqslv)

bayesprev_map <- function(k, n, a=0.05, b=1) {
	# Bayesian maximum a posteriori estimate of population prevalence gamma
	# under a uniform prior
	# 
	# Args:
	#  k: number of participants/tests significant out of 
	#  n: total number of participants/tests
	#  a: alpha value of within-participant test (default=0.05)
	#  b: sensitivity/beta of within-participant test (default=1)
	   
	gm <- (k/n -a)/(b-a)
	if(gm <0) gm <- 0
	if(gm>1) gm <- 1
	return(gm)
} 

bayesprev_posterior <- function(x, k, n, a=0.05, b=1) {
	# Bayesian posterior of population prevalence gamma under a uniform prior
	#
    # Args:
	# x : values of gamma at which to evaluate the posterior density
    # k : number of participants significant out of 
    # n : total number of participants
    # a : alpha value of within-participant test (default=0.05)
    # b : sensitivity/beta of within-participant test (default=1)
    
    m1 <-  k + 1
	m2 <- n - k + 1
    theta <- a + (b-a)*x
    post <- (b -a)*dbeta(theta,m1, m2)
    post <- post/(pbeta(b, m1, m2) - pbeta(a, m1, m2))
    return(post)
}


bayesprev_bound <- function(p, k, n, a=0.05, b=1) {
	# Bayesian lower bound of population prevalence gamma under a uniform prior
	#
    # Args:
	#  p : density the lower bound should bound (e.g. 0.95)
    #  k : number of participants significant out of 
    #  n : total number of participants
    #  a : alpha value of within-participant test (default=0.05)
    #  b : sensitivity/beta of within-participant test (default=1)

	m1 <-  k + 1
	m2 <- n - k + 1
	th_c <- qbeta( p*pbeta(a, m1, m2) + (1-p)*pbeta(b, m1, m2), m1, m2 )
	g_c <- (th_c -a)/(b-a)
	return(g_c)
}


bayesprev_hpdi <- function(p, k, n, a=0.05, b=1) {
	# Bayesian highest posterior density interval of population prevalence gamma
    # under a uniform prior
	#
	# Args:
    #  p : HPDI to return (e.g. 0.95 for 95%)
    #  k : number of participants significant out of 
    #  n : total number of participants
    #  a : alpha value of within-participant test (default=0.05)
    #  b : sensitivity/beta of within-participant test (default=1)

	m1 <- k+1
	m2 <- n-k+1
	
	if(m1 ==1) {
	  endpts <- c(a, qbeta( (1 -p)*pbeta(a, m1, m2) + p*pbeta(b, m1, m2), m1, m2 ) ) 
	  return((endpts -a)/(b-a))
	}
	 
	if(m2 ==1) {
	  endpts <- c( qbeta( p*pbeta(a, m1, m2) + (1- p)* pbeta(b, m1, m2), m1, m2 ) , b) 
	  return( (endpts-a)/(b-a))
	}
	
	if(k<= n*a) {
	  endpts <- c(a, qbeta( (1 -p)*pbeta(a, m1, m2) + p*pbeta(b, m1, m2), m1, m2 ) ) 
	  return((endpts -a)/(b-a))
	}
	 
	if(k>= n*b) {
	  endpts <- c( qbeta( p*pbeta(a, m1, m2) + (1- p)* pbeta(b, m1, m2), m1, m2 ) , b) 
	  return( (endpts-a)/(b-a))
	}


	  g <- function(x, m1, m2, a, b, p ) {
	  y <- numeric(2)
	  y[1] <-  pbeta(x[2], m1, m2) - pbeta(x[1], m1, m2) - p*(pbeta(b, m1, m2) - pbeta(a, m1, m2))
	  y[2] <- log(dbeta(x[2], m1, m2)) - log(dbeta(x[1], m1, m2))
	  return(y)
	}
	   
	x_init <- numeric(2)
	   
	p1 <- (1-p)/2 
	p2 <- (1 +p)/2
	   
	x_init[1] <- qbeta( (1 -p1)*pbeta(a, m1, m2) + p1* pbeta(b, m1, m2), m1, m2 )
	x_init[2] <- qbeta( (1 -p2)*pbeta(a, m1, m2) + p2* pbeta(b, m1, m2), m1, m2 )
	   
	opt <- nleqslv(x_init, g, method ="Newton", control=list(maxit=1000), m1=m1, m2=m2, a=a, b=b, p=p)

	if (opt$termcd ==1)  print("convergence achieved") 
	if (opt$termcd != 1)  print("failed to converge") 
	
	temp <- opt$x
	if (temp[1] <a) {
		temp[1] <- a
		temp[2] <- qbeta( (1 -p)*pbeta(a, m1, m2) + p* pbeta(b, m1, m2), m1, m2 )
	}
	if (temp[2] > b) {
		temp[1] <- qbeta( p*pbeta(a, m1, m2) + (1-p)* pbeta(b, m1, m2), m1, m2 )
		temp[2] <- b
	}
	endpts <- (temp -a)/(b-a)
	return(endpts) 
}


# Functions for  Bayesian inference of difference in prevalence

bayes_prev_diff_between <- function(k1, n1, k2, n2, p, a = 0.05, b = 1, Nsamp){

# Bayesian posterior inference for the difference in prevalence 
# when the same test is applied to two different groups

# Inputs:
#  k1 : number of participants significant in group 1, out of n1
#  n1 : total number of participants in group 1
#  k2 : number of participants significant in group 2, out of n2
#  n2 : total number of participants in group 2
#  p  : coverage for highest-posterior density interval (in [0 1])
#  a  : alpha value of within-participant test (default=0.05)
#  b  : sensitivity/beta of within-participant test (default=1)
#  Nsamp : number of samples from the posterior

#  Outputs:
#  map    : maximum a posteriori estimate of the difference in prevalence:
#                gamma_1 - gamma_2
#  hpdi   : highest-posterior density interval with coverage p 
#  probGT: estimated posterior probability that the prevalence is higher in group 1
#  logoddsGT: estimated log odds in favour of the hypothesis that the prevalence 
#                         is higher in group 1
#  gpost : ggplot object of posterior distribution of the prevalence difference

#  Load the necessary libraries
  
  library(HDInterval)
  library(ggplot2)
  
# Parameters for uniform priors

  r1 =1 ; s1 =1 ; r2 =1 ; s2 = 1

# Parameters for Beta posteriors

  m11 = k1 + r1 ; m12 = n1 - k1 + s1
  m21 = k2 + r2 ; m22 = n2 - k2 + s2
  
# Generate truncated Beta values for theta_1

  vec1 = runif(Nsamp, pbeta(a, m11, m12), pbeta(b, m11, m12))
  th1 = qbeta(vec1, m11, m12)

# Generate truncated Beta values for theta_2

  vec2 = runif(Nsamp, pbeta(a, m21, m22), pbeta(b, m21, m22))
  th2 = qbeta(vec2, m21, m22)

# Compute vector of estimates of prevalence difference

  delta = (th1 - th2)/ (b-a) 
  
# Estimate the posterior probability, and logodds, that the prevalence is higher for group 1.
# Laplace's rule of succession used to avoid estimates of 0 or 1

  
  probGT = (sum(delta > 0)+1)/(Nsamp +2)
  
  logoddsGT = log(probGT / (1 - probGT))

# Compute the HPD interval, coverage = p

  hpdi = hdi(delta, credMass = p)

# Approximate MAP estimate

  dens = density(delta, bw = "SJ")
  map = dens$x[dens$y == max(dens$y)]

# Produce a plot of the posterior density

fill <- "gold1"
line <- "goldenrod2"
dd =data.frame(delta)

gpost <- ggplot(dd, aes(x=delta)) +
        geom_density(fill = fill, colour = line, bw ="SJ") +
        scale_x_continuous(name = "Prevalence difference",
                           limits=c(-1, 1)) +
        scale_y_continuous(name = "Posterior density") +
        ggtitle("Posterior density of difference in Prevalence")

return(list(map = map, hpdi = hpdi, probGT = probGT, logoddsGT = logoddsGT, gpost = gpost))
}


 bayesprev_diff_within <- function(k11, k10, k01, n, p, a = 0.05, b = 1, Nsamp){
 
# Bayesian posterior inference for the difference in prevalence 
# when two different tests are applied to the same group


#Inputs:
#  k11 : number of participants significant in both tests
#  k01 : number of participants significant in test 2 and not test 1
#  k10 : number of participants significant in test 1 and not test 2
#  n   : total number of participants
#  p   : coverage for highest-posterior density interval (in [0 1])
#  a   : alpha value of within-participant test (default=0.05)
#  b   : sensitivity/beta of within-participant test (default=1)
#  Nsamp : number of samples from the posterior

# Outputs:
# map    : maximum a posteriori estimate of the difference in prevalence:
#                 gamma_1 - gamma_2
# hpdi   : highest-posterior density interval with coverage p 
#  probGT: estimated posterior probability that the prevalence is higher on test 1
#  logoddsGT: estimated log odds in favour of the hypothesis that the prevalence 
#                         is higher on test 1
# gpost : ggplot object of the posterior distribution

# Load the necessary libraries

  library(HDInterval)
  library(ggplot2)
  
# Define the parameters  for the Dirichlet prior distribution -- 
# here giving a uniform distribution

  r11 =1; r10 = 1; r01 = 1; r00 = 1

# Compute the parameters for the posterior Dirichlet distribution

  k00 = n - k11 - k10 - k01

  m11 = k11 + r11 ; m10 = k10 + r10 
  m01 = k01 + r01 ; m00 = k00 + r00
  

# Compute samples from the truncated Dirichlet posterior

  z11 = runif(Nsamp, pbeta(0, m11, m10 + m01 + m00), pbeta(b, m11, m10 + m01 + m00))
  u11 = qbeta(z11, m11, m10 + m01 + m00); th11 = u11

  lo = pmax( (a - th11)/(1-th11), 0); hi = (b - th11)/(1-th11)
  z10 = runif(Nsamp, pbeta(lo, m10,  m01 + m00), pbeta(hi, m10, m01 + m00))
  u10 = qbeta(z10, m10,  m01 + m00); th10 = (1 - th11)*u10


  lo = pmax( (a - th11)/(1-th11 - th10), 0); hi = pmin( (b - th11)/(1-th11- th10), 1 )
  z01 = runif(Nsamp, pbeta(lo, m01,   m00), pbeta(hi, m01,  m00))
  u01 = qbeta(z01, m01,   m00); th01 = (1 -th11 - th10)*u01

  th00 = 1 - th11 -th10 - th01

  theta = cbind(th11, th10, th01, th00)

# Compute vector of estimates of prevalence difference

  delta = (th10 - th01)/(b-a)  
  
# Estimate the posterior probability, and logodds, that the prevalence is higher on test 1.
# Laplace's rule of succession used to avoid estimates of 0 or 1
  
  probGT = (sum(delta > 0) +1)/(Nsamp +2)
  
  logoddsGT = log(probGT / (1 - probGT))


# Compute the 96% HPDI

  hpdi = hdi(delta, credMass = p)

  dens = density(delta, bw = "SJ")

# Approximate MAP estimate

  dens = density(delta, bw = "SJ")
  map = dens$x[dens$y == max(dens$y)]

# Produce a plot of the posterior density

  fill <- "gold1"
  line <- "goldenrod2"

  dd = data.frame(delta)

  gpost <- ggplot(dd, aes(x=delta)) +
        geom_density(fill = fill, colour = line, bw ="SJ") +
        scale_x_continuous(name = "Prevalence difference",
                           limits=c(-1, 1)) +
        scale_y_continuous(name = "Posterior density") +
        ggtitle("Posterior density of difference in Prevalence")

  return(list(map = map, hpdi = hpdi, probGT = probGT, logoddsGT = logoddsGT, gpost = gpost))
}

sim_binom <- function(g1, n1, g2, n2, a = 0.05, b = 1){

# Simulation of random binomial data for each of two groups 
# who undertake the same test

th1 = a + (b -a) * g1
th2 = a + (b -a) * g2

k1 = rbinom(1, n1, th1)
k2 = rbinom(1, n2, th2)

return(list(k1 = k1, k2 = k2))
}


sim_multinom <- function(g1, g2, r12, n, a = 0.05, b = 1){

#  Simulation of random multinomial data for the case of 
#  prevalence difference between tests, with a single group

# Input:

#  g1 : prevalence of the effect associated with test 1
#  g2 : prevalence of the effect associated with test 2
#  r12 : correlation between the presence of the effect on test 1 and
#          the presence of the effect on test 2.
#  n   : total number of participants
#  a   : alpha value of within-participant test (default=0.05)
#  b   : sensitivity/beta of within-participant test (default=1)

# Output:

# mdat : a vector containing the numbers of participants who provide
#             a significant result on both tests, on only the first test, 
#             on only the second test, on neither test, respectively
           
           
# Compute the prevalence of the effect on both tests

  g11 = g1 * g2 + r12 * sqrt( g1 *(1-g1) *g2 *(1-g2))
  g10 = g1 - g11
  g01 = g2 - g11
  g00 = 1 - g11 - g10 - g01
  
# Compute the parameters of the Multinomial distribution

  th11 = b^2 * g11 + a * b * g10 + a * b * g01 + a^2 * g00
  th10 = a + (b-a)* g1 - th11
  th01 = a + (b-a)* g2 - th11
  th00= 1 - th11 - th10 -th01
  
  theta =c(th11, th10, th01, th00)
  
  mdat = as.vector(rmultinom(1, n, theta))
  
  return(list(k11 = mdat[1], k10 = mdat[2], k01 = mdat[3], k00 = mdat[4]))
}


