function lo = bayesprev_logodds(k, n, x, a, b)
% Posterior log-odds in favor of the population prevalence gamma being
% greater than x
%
% k : number of participants significant out of 
% n : total number of participants
% x : log-odds threshold (default=0.5)
% a : alpha value of within-participant test (default=0.05)
% b : sensitivity/beta of within-participant test (default=1)

if nargin<=4
    b = 1;
end
if nargin<=3
    a = 0.05;
end
if nargin<=2
    x = 0.5;
end

% gamma prior = Beta(r,s)
r = 1;
s = 1;

the = a + (b-a)*x;

% d = makedist('Beta','a',k+r,'b',n-k+s);
% d.truncate(a,b);
% p = 1-d.cdf(the);

b1 = k+r;
b2 = n-k+s;
p = (betacdf(b,b1,b2)-betacdf(the,b1,b2))./(betacdf(b,b1,b2)-betacdf(a,b1,b2));
lo = log(p/(1-p));
