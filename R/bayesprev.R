map <- function(k, n, a=0.05, b=1) {
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

posterior <- function(x, k, n, a=0.05, b=1) {
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


bound <- function(p, k, n, a=0.05, b=1) {
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
	th_c <- qbeta( (1 -p)*pbeta(a, m1, m2) + p* pbeta(b, m1, m2), m1, m2 )
	g_c <- (th_c -a)/(b-a)
	return(g_c)
}


hdpi <- function(p, k, n, a=0.05, b=1) {
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
	  endpts <- c(a, qbeta( (1 -prob)*pbeta(a, m1, m2) + prob*pbeta(b, m1, m2), m1, m2 ) ) 
	  return((endpts -a)/(b-a))
	}
	 
	if(m2 ==1) {
	  endpts <- c( qbeta( prob*pbeta(a, m1, m2) + (1- prob)* pbeta(b, m1, m2), m1, m2 ) , b) 
	  return( (endpts-a)/(b-a))
	}

	  g <- function(x, m1, m2, a, b, prob ) {
	  y <- numeric(2)
	  y[1] <-  pbeta(x[2], m1, m2) - pbeta(x[1], m1, m2) - prob*(pbeta(b, m1, m2) - pbeta(a, m1, m2))
	  y[2] <- log (dbeta(x[2], m1, m2)) - log(dbeta(x[1], m1, m2))
	  return(y)
	}
	   
	x_init <- numeric(2)
	   
	p1 <- (1-prob)/2 
	p2 <- (1 +prob)/2
	   
	x_init[1] <- qbeta( (1 -p1)*pbeta(a, m1, m2) + p1* pbeta(b, m1, m2), m1, m2 )
	x_init[2] <- qbeta( (1 -p2)*pbeta(a, m1, m2) + p2* pbeta(b, m1, m2), m1, m2 )
	   
	opt <- nleqslv(x_init, g, control=list(maxit=150), m1=m1, m2=m2, a=a, b=b, p=p)

	if (opt$termcd ==1)  print("convergence achieved") 
	if (opt$termcd != 1)  print("failed to converge") 
	
	temp <- opt$x
	if (temp[1] <a) {
		temp[1] <- a
		temp[2] <- qbeta( (1 -prob)*pbeta(a, m1, m2) + prob* pbeta(b, m1, m2), m1, m2 )
	}
	if (temp[2] > b) {
		temp[1] <- qbeta( prob*pbeta(a, m1, m2) + (1-prob)* pbeta(b, m1, m2), m1, m2 )
		temp[2] <- b
	}
	endpts <- (temp -a)/(b-a)
	return(endpts) 
}
