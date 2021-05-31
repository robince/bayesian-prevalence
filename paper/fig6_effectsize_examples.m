% Ince, Paton, Kay and Schyns
% "Bayesian inference of population prevalence"
% biorxiv: https://doi.org/10.1101/2020.07.08.191106

% Figure 6: Examples of different effect size prevalence curves. 
% A,B,C: N=50 participant effects are drawn from a single Gaussian 
% distribution with μ_pop=5,10,7 and σ_b=2,2,0 respectively. All three 
% simulated data sets have the same posterior population prevalence 
% (left hand panels, dot = MAP, thick line = 50% HPDI, 
% thin line = 96% HPDI), but the prevalence effect size curves show 
% different profiles (right hand panels, black line = MAP, grey shaded 
%region = 96% HPDI). D: N=50 participant effects are drawn from a single 
%Gaussian with μ_pop=1, σ_b=2. E,F: 25 participants are simulated with no 
% effect. The other 25 participants are drawn from a Gaussian distribution 
% with E: μ_pop=5,σ_b=3 and F: μ_pop=10,σ_b=1. D,E,F: All three simulated 
% data sets have the same posterior population prevalence (left hand 
% panels), but the prevalence effect size curves show different profiles 
% (right hand panels). k=500 for all simulations. Orange lines show the 
% effect size corresponding to a one-sided α=0.05 within-participant 
%threshold. 


% close all
figure

Ngrp1 = 50;
Ngrp2 = 50;

sigma_w = 10;
Nsamp = 500;

p = 0.05;

load fig6seed
rng(s);

% group 1 
% all members show an effect
% narrow between participant variance
% medium effect size
ax = subplot(3,3,[2 3]);
g1_sig_b = 2;
g1_mu_effect = 5;
rawdat = generate_data(g1_mu_effect, g1_sig_b, sigma_w, Nsamp, Ngrp1);
datA = do_stats(rawdat, Nsamp, Ngrp1);
do_plot(datA,Nsamp,g1_mu_effect,sigma_w,1);
k1 = sum(datA.indt>tinv(1-p/2, 49));
title(sprintf('group 1, k=%d',k1));
% ylabel('Prevalence')


% group 2
% only 25% of members show an effect
% narrow between participant variance for this effect
% strong effect size
ax = subplot(3,3,[5 6]);
g2_sig_b = 2;
g2_mu_effect1 = 2;
g2_mu_effect2 = 10;
g2_prev = 0;
g2_Neff = round(Ngrp2*g2_prev);
g2_Nnoeff = Ngrp2 - g2_Neff;
rawdat1 = generate_data(g2_mu_effect1, g2_sig_b, sigma_w, Nsamp, g2_Neff);
rawdat2 = generate_data(g2_mu_effect2, g2_sig_b, sigma_w, Nsamp, g2_Nnoeff);
rawdat = cat(2,rawdat1,rawdat2);
datB = do_stats(rawdat, Nsamp, Ngrp1);
do_plot(datB,Nsamp,g2_mu_effect2,sigma_w,1);
k2 = sum(datB.indt>tinv(1-p/2, 49));
title(sprintf('group 2, k=%d',k2));
% ylabel('Prevalence')

% group 3
ax = subplot(3,3,[8 9]);
g3_sig_b = 0;
g3_mu_effect = 7;
rawdat = generate_data(g3_mu_effect, g3_sig_b, sigma_w, Nsamp, Ngrp1);
datC = do_stats(rawdat, Nsamp, Ngrp1);
do_plot(datC,Nsamp,g1_mu_effect,sigma_w,1);
k3 = sum(datC.indt>tinv(1-p/2, 49));
title(sprintf('group 3, k=%d',k3));
xlabel(sprintf('T(%d)',Nsamp-1));


% %%
% figure
ax = [];
ax(1) = subplot(3,3,1);
prev_plot(k1,Ngrp1);
ylabel('Prevalence')

ax(2) = subplot(3,3,4);
prev_plot(k2,Ngrp1);
ylabel('Prevalence')

ax(3) = subplot(3,3,7);
prev_plot(k3,Ngrp1);
ylabel('Prevalence')

set(ax,'XTick',1)
set(ax,'XTickLabel',[])
set(ax,'ylim',[0.8 1])


%%
figure


% group 1 
% 50% prevalence single dist
ax = subplot(3,3,[2 3]);
g1_sig_b = 2;
g1_mu_effect = 1;
rawdat = generate_data(g1_mu_effect, g1_sig_b, sigma_w, Nsamp, Ngrp1);
datA = do_stats(rawdat, Nsamp, Ngrp1);
do_plot(datA,Nsamp,g1_mu_effect,sigma_w,1);
k1 = sum(datA.indt>tinv(1-p/2, 49));
title(sprintf('group 1, k=%d',k1));
% ylabel('Prevalence')


% group 2
ax = subplot(3,3,[5 6]);
g2_sig_b = 3;
g2_mu_effect1 = 0;
g2_mu_effect2 = 5;
g2_prev = 0.5;
g2_Neff = round(Ngrp2*g2_prev);
g2_Nnoeff = Ngrp2 - g2_Neff;
rawdat1 = generate_data(g2_mu_effect1, 0, sigma_w, Nsamp, g2_Nnoeff);
rawdat2 = generate_data(g2_mu_effect2, g2_sig_b, sigma_w, Nsamp, g2_Neff);
rawdat = cat(2,rawdat1,rawdat2);
datB = do_stats(rawdat, Nsamp, Ngrp1);
do_plot(datB,Nsamp,g2_mu_effect2,sigma_w,1);
k2 = sum(datB.indt>tinv(1-p/2, 49));
title(sprintf('group 2, k=%d',k2));
% ylabel('Prevalence')

% group 3
ax = subplot(3,3,[8 9]);
g2_sig_b = 1;
g2_mu_effect1 = 0;
g2_mu_effect2 = 10;
g2_prev = 0.5;
g2_Neff = round(Ngrp2*g2_prev);
g2_Nnoeff = Ngrp2 - g2_Neff;
rawdat1 = generate_data(g2_mu_effect1, 0, sigma_w, Nsamp, g2_Nnoeff);
rawdat2 = generate_data(g2_mu_effect2, g2_sig_b, sigma_w, Nsamp, g2_Neff);
rawdat = cat(2,rawdat1,rawdat2);
datC = do_stats(rawdat, Nsamp, Ngrp1);
do_plot(datC,Nsamp,g2_mu_effect2,sigma_w,1);
k3 = sum(datC.indt>tinv(1-p/2, 49));
title(sprintf('group 3, k=%d',k2));
xlabel(sprintf('T(%d)',Nsamp-1));


% %%
% figure
ax = [];
ax(1) = subplot(3,3,1);
prev_plot(k1,Ngrp1);
ylabel('Prevalence')

ax(2) = subplot(3,3,4);
prev_plot(k2,Ngrp1);
ylabel('Prevalence')

ax(3) = subplot(3,3,7);
prev_plot(k3,Ngrp1);
ylabel('Prevalence')

set(ax,'XTick',1)
set(ax,'XTickLabel',[])
set(ax,'ylim',[0 1])

%%


function prev_plot(k,Nsub)
i=1;
oil = 3;
iil = 8;
hy = 0.3;
a = 0.05;
b = 1;
co = get(gca,'ColorOrder');
lw = 2;
hold on

xmap = bayesprev_map(k, Nsub, a, b);
pmap = bayesprev_posterior(xmap, k, Nsub, a, b);        %#ok<NASGU>

yp = 1;
plot(yp, xmap, '.','MarkerSize',20,'Color',co(i,:));

h = bayesprev_hpdi(0.96,k, Nsub, a, b);
plot([yp yp],[h(1) h(2)],'Color',[co(i,:) 0.3],'LineWidth',oil)

h = bayesprev_hpdi(0.5,k,Nsub, a, b);
plot([yp yp],[h(1) h(2)],'Color',[co(i,:) 0.3],'LineWidth',iil)


end

function ax = do_plot(dat,Nsamp,mu_g,sigma_w,pi)
[es pmap hpdi] = prev_curve_onesided(dat,1);
posbar = hpdi(2,:) - pmap;
negbar = pmap - hpdi(1,:);
shadedErrorBar(es,pmap,cat(1,posbar,negbar))
% xline(mu_g./(sigma_w./sqrt(Nsamp)),'b')
p = 0.05;
xline(tinv(1-p/2, Nsamp-1),'r');
xline(tinv(p/2,Nsamp-1),'r');
ylim([0 1])
xlim([0 30])
% title('Right-tailed prevalence')


% ax(2) = subplot(n,m,pi+m);
% [es pmap hpdi] = prev_curve_onesided(dat,-1);
% posbar = hpdi(2,:) - pmap;
% negbar = pmap - hpdi(1,:);
% shadedErrorBar(es,pmap,cat(1,posbar,negbar))
% title('Left-tailed prevalence')
% xline(mu_g./(sigma_w./sqrt(Nsamp)),'b')
% p = 0.05;
% xline(tinv(1-p/2, Nsamp-1),'r');
% xline(tinv(p/2,Nsamp-1),'r');
% linkaxes(ax,'xy')
% ylim([0 0.8])
end

% generate data from heirachical normal model
function dat = generate_data(mu_g, sigma_b, sigma_w, Nsamp, Nsub)
% mu_g - ground truth group mean
% sigma_b - between participant standard deviation
% sigma_w - within participant standard deviation
% Nsamp - number of trials per participant
% Nsub - number of participants

% generate individual subject means from population normal distribution
submeanstrue = normrnd(mu_g, sigma_b, [Nsub 1]);
dat = zeros(Nsamp, Nsub);
for si=1:Nsub
    % generate within-participant data
    dat(:,si) = normrnd(submeanstrue(si), sigma_w, [Nsamp 1]);
end
end

function dat = do_stats(rawdat, Nsamp, Nsub)
% generate data from heirachical normal model
dat = [];
dat.Nsub = Nsub;
dat.Nsamp = Nsamp;
dat.indsig = false(1,Nsub);
dat.indt = zeros(1,Nsub);
for si=1:Nsub
    % within-participant t-test significance
    [dat.indsig(si) p ci stats] = ttest(rawdat(:,si));
    % within-participant t-score
    dat.indt(si) = stats.tstat;
end
end