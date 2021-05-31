%
% Ince, Paton, Kay and Schyns
% "Bayesian inference of population prevalence"
% biorxiv: https://doi.org/10.1101/2020.07.08.191106
%
% Figure 4: Example where between-group prevalence diverges from two-sample 
% t-test. 
% We simulate standard hierarchical Gaussian data for two groups of 20 
% participants, T=100,  σ_w=10, N=20 per group. A: Group 1 participants are 
% drawn from a single population Gaussian distribution with μ=4,σ_b=1. 
% Group 2 participants are drawn from two Gaussian distributions. 75% of 
% participants are drawn from N(0,0.01) and 25% of participants are drawn 
% from N(16,0.5). Dashed line shows with the p=0.05 within-participant 
% threshold (one-sample t-test). The means of these two groups are not 
% significantly different (B), but they have very different prevalence 
% posteriors (C). The posterior distribution for the difference in 
% prevalence shows the higher prevalence in Group 1: 
% 0.61 [0.36 0.85] (MAP [96% HPDI]) (D).  

Ngrp1 = 20;
Ngrp2 = 20;

sig_w = 10;
Nsamp = 100;

% group 1 
% all members show an effect
% narrow between participant variance
% medium effect size
g1_sig_b = 1;
g1_mu_effect = 4;
g1dat = generate_data(g1_mu_effect, g1_sig_b, sig_w, Nsamp, Ngrp1);

% group 2
% only 25% of members show an effect
% narrow between participant variance for this effect
% strong effect size
g2_sig_b = 0.5;
g2_mu_effect = 16;
g2_prev = 0.25;
g2_Neff = round(Ngrp2*g2_prev);
g2_Nnoeff = Ngrp2 - g2_Neff;
g2dat_effect = generate_data(g2_mu_effect, g2_sig_b, sig_w, Nsamp, g2_Neff);
g2dat_noeffect = generate_data(0, 0.01, sig_w, Nsamp, g2_Nnoeff);
g2dat = cat(2,g2dat_effect,g2dat_noeffect);

% between group t-test
[tsig grp_p ci stats] = ttest2(mean(g1dat),mean(g2dat));
grp_t = stats.tstat;

% within participant tests
g1_indsig = false(1,Ngrp1);
g1_indt = zeros(1,Ngrp1);
for si=1:Ngrp1
    % within-participant t-test significance
    [g1_indsig(si) p ci stats] = ttest(g1dat(:,si));
    % within-participant t-score
    g1_indt(si) = stats.tstat;
end
g2_indsig = false(1,Ngrp2);
g2_indt = zeros(1,Ngrp2);
for si=1:Ngrp2
    % within-participant t-test significance
    [g2_indsig(si) p ci stats] = ttest(g2dat(:,si));
    % within-participant t-score
    g2_indt(si) = stats.tstat;
end

[map, px, p, hpdi, probGT, loGT, samples] = bayesprev_diff_between(sum(g1_indsig), Ngrp1, sum(g2_indsig), Ngrp2, 0.96);

% plots
figure
mg1 = mean(g1dat);
mg2 = mean(g2dat);

subplot(1,3,1)
violins = violinplot([mg1' mg2'],{'Group 1','Group 2'});
% title('Participant means in each group')
c1 = violins(1).ViolinColor;
c2 = violins(2).ViolinColor;
hold on
yline(tinv(1-0.05,99),'k--')

subplot(1,3,2)
hold on
h = bar(1,mean(mg1),'k');
set(h,'facecolor',c1)
set(h,'facealpha',violins(1).ViolinAlpha)
h = bar(2,mean(mg2),'k');
set(h,'facecolor',c2);
set(h,'facealpha',violins(1).ViolinAlpha)

er = errorbar(1 , mean(mg1), std(mg1)./sqrt(Ngrp1) );
er.LineStyle = 'none';
er.Color = 'k';
er.LineWidth = 2;
set(er.Bar, 'ColorType', 'truecoloralpha', 'ColorData', [er.Line.ColorData(1:3); 255*0.3])

er = errorbar(2, mean(mg2), std(mg2)/sqrt(Ngrp2));
er.LineStyle = 'none';
er.Color = 'k';
er.LineWidth = 2;
set(er.Bar, 'ColorType', 'truecoloralpha', 'ColorData', [er.Line.ColorData(1:3); 255*0.3])
ylabel('Mean')
% title('T-test not significant p=0.92')
set(gca,'XTick', [1 2])
set(gca,'XTickLabels',{'Group 1', 'Group 2'})


subplot(2,3,3)
box off
xlabel('Prevalence Proportion')
ylabel('Posterior Density')
hold on
a = 0.05;
b = 1;

oil = 2;
iil = 8;
hy = 0.3;
lw = 2;
alf = violins(1).ViolinAlpha;
% co = get(gca,'ColorOrder');
% co = [c1;c2];
co = cat(1,[c1 alf], [c2 alf]);
x = linspace(0,1,100);

k = sum(g1_indsig);
N = Ngrp1;
i=1;

lh(1) = plot(x, bayesprev_posterior(x, k, N, a, b),'Color',co(i,:),'LineWidth',lw);
xmap = bayesprev_map(k, N, a, b);
pmap = bayesprev_posterior(xmap, k, N, a, b);        %#ok<NASGU>
h = bayesprev_hpdi(0.96,k, N, a, b);

yp = 0.5;
yp = 0.25;
plot(xmap, yp,'.','MarkerSize',20,'Color',co(i,:));

plot([h(1) h(2)],[yp yp],'Color',co(i,:),'LineWidth',oil)
h = bayesprev_hpdi(0.5,k,N, a, b);
plot([h(1) h(2)],[yp yp],'Color',co(i,:),'LineWidth',iil)

k = sum(g2_indsig);
N = Ngrp2;
i=2;

lh(1) = plot(x, bayesprev_posterior(x, k, N, a, b),'Color',co(i,:),'LineWidth',lw);

xmap = bayesprev_map(k, N, a, b);
pmap = bayesprev_posterior(xmap, k, N, a, b);        %#ok<NASGU>
h = bayesprev_hpdi(0.96,k, N, a, b);

% yp = pmap;
yp = 0.5;
yp = 0.25;
plot(xmap, yp,'.','MarkerSize',20,'Color',co(i,:));

plot([h(1) h(2)],[yp yp],'Color',co(i,:),'LineWidth',oil)
h = bayesprev_hpdi(0.5,k,N, a, b);
plot([h(1) h(2)],[yp yp],'Color',co(i,:),'LineWidth',iil)

subplot(2,3,6)
hold on
plot(px,p,'color',[0 0 0 0.3],'LineWidth',lw)
yp = 0.25;
plot(map, yp,'.','MarkerSize',20,'Color',[0 0 0 0.3]);

plot([hpdi(1) hpdi(2)],[yp yp],'Color',[0 0 0 0.3],'LineWidth',iil)
box off
xlabel('Prevalence Difference \gamma_1 - \gamma_2')
ylabel('Posterior Density')


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