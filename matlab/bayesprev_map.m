function gamma = bayesprev_map(k, n, a, b)
% Bayesian maximum a posteriori estimate of population prevalence gamma
% under a uniform prior
%
% k : number of participants significant out of 
% n : total number of participants
% a : alpha value of within-participant test (default=0.05)
% b : sensitivity/beta of within-participant test (default=1)

if nargin<=3
    b = 1;
end
if nargin<=2
    a = 0.05;
end

% gamma prior = Beta(r,s)
r = 1;
s = 1;

theta = (k+r-1)./(n+r+s-2);
if theta<=a
    gamma = 0;
elseif theta>=b
    gamma = 1;
else
    gamma = (theta-a)/(b-a);
end

