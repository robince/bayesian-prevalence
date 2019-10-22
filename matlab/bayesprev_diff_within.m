function [map, post_x, post, hpi] = bayesprev_diff_within(k11, k01, k10, n, p, a, b, Nsamp)
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
the11d = the11d.truncate(0,b);

the10d = makedist('Beta','a',m10,'b',m01+m00);
the01d = makedist('Beta','a',m01,'b',m00);

theta_samp = zeros(4,Nsamp);

% sample from truncated dirichlet
for si=1:Nsamp
    the11 = the11d.random(1,1);
    theta_samp(1,si) = the11;
    
    the10d_trunc = the10d.truncate( max([ (a-the11)./(1-the11) 0 ]), (b-the11)./(1-the11) );
    u10 = the10d_trunc.random(1,1);
    the10 = (1-the11).*u10;
    theta_samp(2,si) = the10;
    
    the01d_trunc = the01d.truncate( max([ (a-the11)./(1-the11-the10) 0 ]), min([ (b-the11)./(1-the11-the10) 1]) );
    the01 = (1-the11-the10).*the01d_trunc.random(1,1);
    theta_samp(3,si) = the01;
end

g_diff_samples = (theta_samp(2,:) - theta_samp(3,:))./(b-a);

% [min(g_diff_samples) max(g_diff_samples)]

x = linspace(-1,1,200); 
g_diff_post = ksdensity(g_diff_samples,x);

[~, idx] = max(g_diff_post);
g_diff_map = x(idx);

map = g_diff_map;
post = g_diff_post;
post_x = x;
if nargin>=5
    hpi = hpdi(g_diff_samples',100*p);
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
[foo,hpdpi]=min(cis);
if m==1
  hpdi=[pt1(hpdpi); pt2(hpdpi)];
else
  hpdpi=sub2ind(size(pt1),hpdpi,1:m);
  hpdi=[pt1(hpdpi); pt2(hpdpi)];
end
