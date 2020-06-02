function [map, post_x, post_p, hpi, post_samples] = bayesprev_diff_within(k11, k01, k10, n, p, a, b, Nsamp)
% Bayesian maximum a posteriori estimate of the difference in prevalence 
% when two tests are applied to the same group
%
% k11 : number of participants significant in both tests
% k01 : number of participants significant in test 2 and not test 1
% k10 : number of participants significant in test 1 and not test 2
% n   : total number of participants
% p   : coverage for highest-posterior density interval (in [0 1])
% a   : alpha value of within-participant test (default=0.05)
% b   : sensitivity/beta of within-participant test (default=1)
% Nsamp : number of samples from the posterior
%
% Outputs:
% map    : maximum a posteriori estimate of the difference in prevalence:
%          gamma_1 - gamma_2
% post_x : x-axis for kernel density fit of posterior distribution of the
%          above
% post   : posterior distribution from kernel density fit
% hpdi   : highest-posterior density interval with coverage p 
if nargin<=6
    b = 1;
end
if nargin<=5
    a = 0.05;
end
if nargin<=7
    Nsamp = 10000;
end

% default both tests have the same properties
if length(a)==1
    a = repmat(a,1,2);
end
if length(b)==1
    b = repmat(b,1,2);
end

% gamma priors, Dirichlet parameters
r11 = 1;
r10 = 1;
r01 = 1;
r00 = 1;

% posterior dirichlet parameters
k00 = n - k11 - k01 - k10;

if k00<0
    error("Input test results don't sum to n")
end

m11 = k11 + r11;
m10 = k10 + r10;
m01 = k01 + r01;
m00 = k00 + r00;

the11d = makedist('Beta','a',m11,'b',m01+m10+m00);
the11d = the11d.truncate(0,min(b));

the10d = makedist('Beta','a',m10,'b',m01+m00);
the01d = makedist('Beta','a',m01,'b',m00);

the11 = the11d.random(Nsamp,1);

low = max((a(1)-the11)./(1-the11),0);
high = (b(1)-the11)./(1-the11);
zlow = the10d.cdf(low);
zhigh = the10d.cdf(high);
z10 = zlow + (zhigh-zlow).*rand(Nsamp,1);
u10 = the10d.icdf(z10);
the10 = (1-the11).*u10;

low = max((a(2)-the11)./(1-the11-the10), 0);
high = min((b(2)-the11)./(1-the11-the10), 1);
zlow = the01d.cdf(low);
zhigh = the01d.cdf(high);
z01 = zlow + (zhigh-zlow).*rand(Nsamp,1);
u01 = the01d.icdf(z01);
the01 = (1-the11-the10).*u01;

g_diff_samples = (the11 + the10 - a(1))./(b(1)-a(1)) - (the11 + the01 - a(2))./(b(2)-a(2));

x = linspace(-1,1,200); 
g_diff_post = ksdensity(g_diff_samples,x);

[~, idx] = max(g_diff_post);
g_diff_map = x(idx);

map = g_diff_map;
post_p = g_diff_post;
post_x = x;
post_samples = g_diff_samples;
if nargin>=5
    hpi = hpdi(g_diff_samples,100*p);
end


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
