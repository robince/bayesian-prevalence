% Ince, Kay and Schyns
% "Bayesian inference of population prevalence"
% biorxiv: https://doi.org/10.1101/2020.07.08.191106
%
% Figure 4: One-sided prevalence as a function of effect size.    

Nsub = 50;
sigma_w = 10;
sigma_b = 2;

figure
ax = [];

% s = rng;
load('currfig4seed.mat')
rng(s)
% MODEL A
Nsamp = 20;
mu_g = 0;
sigma_g = sqrt(sigma_b.^2 + ((sigma_w).^2)/Nsamp);
dat = generate_data(mu_g, sigma_b, sigma_w, Nsamp, Nsub);
dat.sigma_g = sigma_g;
ax = do_plot(dat,Nsamp,mu_g,sigma_w,1)

% MODEL B
Nsamp = 500;
mu_g = 0;
sigma_g = sqrt(sigma_b.^2 + ((sigma_w).^2)/Nsamp);
dat = generate_data(mu_g, sigma_b, sigma_w, Nsamp, Nsub);
dat.sigma_g = sigma_g;
ax = do_plot(dat,Nsamp,mu_g,sigma_w,2)

% MODEL C
Nsamp = 20;
mu_g = 1;
sigma_g = sqrt(sigma_b.^2 + ((sigma_w).^2)/Nsamp);
dat = generate_data(mu_g, sigma_b, sigma_w, Nsamp, Nsub);
dat.sigma_g = sigma_g;
ax = do_plot(dat,Nsamp,mu_g,sigma_w,5)

% MODEL D
Nsamp = 500;
mu_g = 1;
sigma_g = sqrt(sigma_b.^2 + ((sigma_w).^2)/Nsamp);
dat = generate_data(mu_g, sigma_b, sigma_w, Nsamp, Nsub);
dat.sigma_g = sigma_g;
ax = do_plot(dat,Nsamp,mu_g,sigma_w,6)

function ax = do_plot(dat,Nsamp,mu_g,sigma_w,pi)
n = 4;
m = 2;
ax(1) = subplot(n,m,pi);
[es pmap hpdi] = prev_curve_onesided(dat,1);
posbar = hpdi(2,:) - pmap;
negbar = pmap - hpdi(1,:);
shadedErrorBar(es,pmap,cat(1,posbar,negbar))
xline(mu_g./(sigma_w./sqrt(Nsamp)),'b')
p = 0.05;
xline(tinv(1-p/2, Nsamp-1),'r');
xline(tinv(p/2,Nsamp-1),'r');
ylim([0 0.8])


title('Right-tailed prevalence')
ax(2) = subplot(n,m,pi+m);
[es pmap hpdi] = prev_curve_onesided(dat,-1);
posbar = hpdi(2,:) - pmap;
negbar = pmap - hpdi(1,:);
shadedErrorBar(es,pmap,cat(1,posbar,negbar))
title('Left-tailed prevalence')
xline(mu_g./(sigma_w./sqrt(Nsamp)),'b')
p = 0.05;
xline(tinv(1-p/2, Nsamp-1),'r');
xline(tinv(p/2,Nsamp-1),'r');
linkaxes(ax,'xy')
ylim([0 0.8])
end
