function hpdi = bayesprev_hpdi(p, k, n, a, b)
% Bayesian highest posterior density interval of population prevalence gamma
% under a uniform prior
%
% p : HPDI to return (e.g. 0.95 for 95%)
% k : number of participants significant out of 
% n : total number of participants
% a : alpha value of within-participant test (default=0.05)
% b : sensitivity/beta of within-participant test (default=1)

if nargin<=5
    b = 1;
end
if nargin<=4
    a = 0.05;
end
if nargin<=3
    inc = 0.01;
end

% gamma prior = Beta(r,s)
r = 1;
s = 1;

td = makedist('Beta','a',k+r,'b',n-k+s);
td = td.truncate(a,b);

if k==0
    x = [a td.icdf(p)];
elseif k==n
    x = [td.icdf(1-p) b];
else
    f = @(x) [td.cdf(x(2))-td.cdf(x(1))-p, td.pdf(x(2))-td.pdf(x(1))];
    opt.Display = 'off';
    opt.FunctionTolerance = 1e-10;
    [x, fval, exitflag, output] = fsolve(f, [td.icdf((1-p)/2) td.icdf((1+p)/2)],opt);
end

% limit to valid theta values
if x(1)<a || x(2)<x(1)
    x = [a td.icdf(p)];
end
if x(2)>b 
    x = [td.icdf(1-p) b];
end
hpdi = (x-a)./(b-a);
