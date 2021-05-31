% Ince, Paton, Kay and Schyns
% "Bayesian inference of population prevalence"
% biorxiv: https://doi.org/10.1101/2020.07.08.191106
%
% Figure 5: One-sided prevalence as a function of effect size. 
% We consider the same simulated systems shown in Figure 1, showing both 
% right-tailed (E_p>E ̂) and left-tailed (E_p<E ̂) prevalence as a function o
% f effect size. Orange lines show the effect size corresponding to the 
% two-sided α=0.05 within-participant test, as used in Figure 1. Dashed 
% lines show the effect size corresponding to the ground truth of the 
% simulation. A,B: μ_pop=0,  C,D: μ_pop=1. A,C: T =20 trials, 
% B,D: T =500 trials. Black line shows MAP, shaded region shows 96% HPDI.   

Nsub = 50;
sigma_w = 10;
sigma_b = 2;

figure
ax = [];

% s = rng;
load('fig5seed.mat')
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
