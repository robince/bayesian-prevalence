function post = bayesprev_posterior(x, k, n, a, b)
% Bayesian posterior of population prevalence gamma
% under a uniform prior
%
% x : values of gamma at which to evaluate the posterior density
% k : number of participants significant out of 
% n : total number of participants
% a : alpha value of within-participant test (default=0.05)
% b : sensitivity/beta of within-participant test (default=1)

if nargin<=4
    b = 1;
end
if nargin<=3
    a = 0.05;
end

if any(x<0) || any(x>1)
    error('prev_posterior: requested value out of range')
end
% gamma prior = Beta(r,s)
r = 1;
s = 1;

theta = a+(b-a).*x;
% d = makedist('Beta','a',k+r,'b',n-k+s);
% d.truncate(a,b);
% post = d.pdf(x);
post = (b-a).*betapdf(theta, k+r, n-k+s);
post = post ./ (betacdf(b, k+r, n-k+s) - betacdf(a, k+r, n-k+s));

