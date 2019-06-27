function g_c = prev_bound(p, k, n, a, b)
% Bayesian lower bound of population prevalence gamma
% under a uniform prior
%
% p : density the lower bound should bound (e.g. 0.95)
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

% gamma prior = Beta(r,s)
r = 1;
s = 1;


% td = makedist('Beta','a',b1,'b',b2);
% td.truncate(a,b);
% the_c = td.icdf(1-p);

% % f = @(the) betacdf(the,b1,b2)-(1-p).*betacdf(b,b1,b2)-p.*betacdf(a,b1,b2);
% opt.Display = 'off';
% start search at mean of posterior
% [the_c, fval, exitflag, output] = fsolve(f, b1./(b1+b2),opt);
% the_c = lsqnonlin(f,b1./(b1+b2),a,b);

% if exitflag<1
%     disp(output.message);
%     error('prev_bound: fsolve failed to solve')
% end

% post = (b-a).*betapdf(theta, b1, b2);
% post = post ./ (betacdf(b, k+r, n-k+s) - betacdf(a, k+r, n-k+s));

b1 = k+r;
b2 = n-k+s;
cdfp = (1-p)*betacdf(b,b1,b2) + p*betacdf(a,b1,b2);
the_c = betainv(cdfp,b1,b2);
g_c = (the_c - a)./(b-a);