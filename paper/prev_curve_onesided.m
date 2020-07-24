function [es pmap hpdi] = prev_curve_onesided(dat,side)
% prevalence of a one-sided effect size threshold for a t-test

Nsamp = dat.Nsamp;
Nsub = dat.Nsub;
Nx = 100;
edat = dat.indt;
esx = linspace(min(edat),max(edat),Nx);
emap = zeros(1,Nx);
eh = zeros(2, Nx);
b = 1;
for xi=1:Nx
%     if xi==24
%         keyboard
%     end
    % number greater than threshold
    if side>0
        k = sum(edat>esx(xi));
        a = 1 - tcdf(esx(xi),Nsamp-1);
    elseif side<0
        k = sum(edat<esx(xi));
        a = tcdf(esx(xi),Nsamp-1);
    end
    if a>=0.8 || abs(b-a) < 2*eps(b)
        emap(xi) = NaN;
        eh(:,xi) = NaN;
        continue
    end
    emap(xi) = bayesprev_map(k,Nsub,a,b);
    try
        eh(:,xi) = bayesprev_hpdi(0.96,k,Nsub,a,b);
    catch
        % a,b too close together, distribution 
        emap(xi) = NaN;
        eh(:,xi) = NaN;
        continue
    end
end

es = esx;
pmap = emap;
hpdi = eh;