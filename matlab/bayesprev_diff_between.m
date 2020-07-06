function [map, post_x, post_p, hpi, probGT, logoddsGT, samples] = bayesprev_diff_between(k1, n1, k2, n2, p, a, b, Nsamp)
% Bayesian maximum a posteriori estimate of the difference in prevalence 
% when the same test is applied to two groups
%
% k1 : number of participants significant in group 1 out of 
% n1 : total number of participants in group 1
% k2 : number of participants significant in group 2 out of 
% n2 : total number of participants in group 2
% p  : coverage for highest-posterior density interval (in [0 1])
% a  : alpha value of within-participant test (default=0.05)
% b  : sensitivity/beta of within-participant test (default=1)
% Nsamp : number of samples from the posterior
%
% Outputs:
% map    : maximum a posteriori estimate of the difference in prevalence:
%          gamma_1 - gamma_2
% post_x : x-axis for kernel density fit of posterior distribution of the
%          above
% post   : posterior distribution from kernel density fit
% hpdi   : highest-posterior density interval with coverage p
% probGT : estimated posterior probability that the prevalence is higher in group 1
% logoddsGT : estimated log odds in favour of the hypothesis that the prevalence is 
%             higher in group 1
% samples : posterior samples
if nargin<=6
    b = 1;
end
if nargin<=5
    a = 0.05;
end
if nargin<=4
    p = 0.96;
end
if nargin<=7
    Nsamp = 10000;
end

% gamma priors = Beta(r,s)
r1 = 1;
s1 = 1;
r2 = 1;
s2 = 1;

the1d = makedist('Beta','a',k1+r1,'b',n1-k1+s1);
the1d = the1d.truncate(a,b);

the2d = makedist('Beta','a',k2+r2,'b',n2-k2+s2);
the2d = the2d.truncate(a,b);

the1samp = the1d.random(Nsamp,1);
the2samp = the2d.random(Nsamp,1);

g_diff_samples = (the1samp - the2samp)./(b-a);

x = linspace(-1,1,200); 
g_diff_post = ksdensity(g_diff_samples,x);

[~, idx] = max(g_diff_post);
g_diff_map = x(idx);

map = g_diff_map;
post_p = g_diff_post;
post_x = x;
samples = g_diff_samples;
hpi = hpdi(g_diff_samples,100*p);
% Estimate the posterior probability, and logodds, that the prevalence is higher for group 1.
% Laplace's rule of succession used to avoid estimates of 0 or 1
probGT = (sum(samples>0)+1)/(Nsamp+2);
logoddsGT = log(probGT / (1-probGT));

function hpdi = hpdi(x, p)
% HPDI - Estimates the Bayesian HPD intervals
%
%   Y = HPDI(X,P) returns a Highest Posterior Density (HPD) interval
%   for each column of X. P must be a scalar. Y is a 2 row matrix
%   where ith column is HPDI for ith column of X.

%   References:
%      [1] Chen, M.-H., Shao, Q.-M., and Ibrahim, J. Q., (2000).
%          Monte Carlo Methods in Bayesian Computation. Springer-Verlag.

% Copyright (C) 2001 Aki Vehtari
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

if nargin < 2
  error('Not enough arguments')
end

m=size(x,2);
pts=linspace(0,100-p,100);
pt1=prctile(x,pts);
pt2=prctile(x,p+pts);
cis=abs(pt2-pt1);
[~,hpdpi]=min(cis);
if m==1
  hpdi=[pt1(hpdpi); pt2(hpdpi)];
else
  hpdpi=sub2ind(size(pt1),hpdpi,1:m);
  hpdi=[pt1(hpdpi); pt2(hpdpi)];
end