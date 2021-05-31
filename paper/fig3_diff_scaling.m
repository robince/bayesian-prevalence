% Ince, Paton, Kay and Schyns
% "Bayesian inference of population prevalence"
% biorxiv: https://doi.org/10.1101/2020.07.08.191106
% 
% Figure 3: Bayesian inference of difference of prevalence. 
% A,B: We consider two independent groups of participants with population 
% prevalence of true positives [γ_1,γ_2] of [25% 25%] (blue), 
% [25% 50%] (red) and [25% 75%] (yellow). We show how A: the Bayesian MAP 
% estimate, and B: 96% HPDI width, of the estimated between-group 
% prevalence difference γ_1-γ_2 scale with the number of participants. 
% C,D: We consider two tests applied to the same sample of participants. 
% Here, each simulation is parameterised by the population prevalence of 
% true positives for the two tests, [γ_1,γ_2], as well as ρ_12, the 
% correlation between the (binary) test results across the population. We 
% show this for [50% 50%] with ρ_12=0.2 (blue), [50% 0%] with ρ_12=0 (red), 
% and [75% 50%] with ρ_12=-0.2 (yellow). We show how C: the Bayesian MAP 
% estimate, and D: 96% HPDI width, of the estimated within-group prevalence 
% difference γ_1-γ_2 scale with the number of participants.

a = 0.05;
b = 1;

figure
hold all
ax = [];
c = get(0, 'DefaultAxesColorOrder');

% 
% Panels A and B
%
% bayes_scale_between.mat from run_scaling_between.m
load bayes_scale_between
ax(1) = subplot(2,2,1);
vi = 1;
gi=1;
dat = squeeze(res(vi,:,gi,:));
shadedErrorBar(Nvals, mean(dat), std(dat),'lineprops',{'-' 'color' c(gi,:)})
gi=2;
dat = squeeze(res(vi,:,gi,:));
shadedErrorBar(Nvals, mean(dat), std(dat),'lineprops',{'-' 'color' c(gi,:)})
gi=3;
dat = squeeze(res(vi,:,gi,:));
shadedErrorBar(Nvals, mean(dat), std(dat),'lineprops',{'-' 'color' c(gi,:)})

% legend({'[25% 25%]' '[25% 50%]' '[25% 75%]'},'location','southeast')
legend({'[0.25 0.25]' '[0.25 0.5]' '[0.25 0.75]'},'location','southeast')
xlabel('Participants (N)')
ylabel('MAP')
title('Between - Bayesian MAP')
axis square
grid off

ax(2) = subplot(2,2,2);
vi = 2;
gi=1;
dat = squeeze(res(vi,:,gi,:));
shadedErrorBar(Nvals, mean(dat), std(dat),'lineprops',{'-' 'color' c(gi,:)})
gi=2;
dat = squeeze(res(vi,:,gi,:));
shadedErrorBar(Nvals, mean(dat), std(dat),'lineprops',{'-' 'color' c(gi,:)})
gi=3;
dat = squeeze(res(vi,:,gi,:));
shadedErrorBar(Nvals, mean(dat), std(dat),'lineprops',{'-' 'color' c(gi,:)})

% legend({'[25% 25%]' '[25% 50%]' '[25% 75%]'},'location','southeast')
legend({'[0.25 0.25]' '[0.25 0.5]' '[0.25 0.75]'},'location','northeast')
xlabel('Participants (N)')
ylabel('HPDI Width')
title('Between - 96% HPDI Width')
axis square
grid off


% 
% Panels C and D
%
% bayes_scale_within.mat from run_scaling_within.m
load bayes_scale_within
gts = [0.5 0.5 0.2; 0 0.5 0; 0.5 0.75 -0.2];
a = 0.05;
b = 1;

ax(3) = subplot(2,2,3);
vi = 1;
gi=1;
dat = squeeze(res(vi,:,gi,:));
shadedErrorBar(Nvals, mean(dat), std(dat),'lineprops',{'-' 'color' c(gi,:)})
gi=2;
dat = squeeze(res(vi,:,gi,:));
shadedErrorBar(Nvals, mean(dat), std(dat),'lineprops',{'-' 'color' c(gi,:)})
gi=3;
dat = squeeze(res(vi,:,gi,:));
shadedErrorBar(Nvals, mean(dat), std(dat),'lineprops',{'-' 'color' c(gi,:)})
% legend({'[0.5 0.5] \rho=0' ...
%     '[0 0.5]\rho=0'...
%     '[0.5 0.75] \rho=0.2'},'location','southeast')
xlabel('Participants (N)')
ylabel('MAP')
title('Within - Bayesian MAP')
axis square
grid off

ax(4) = subplot(2,2,4);
vi = 2;
gi=1;
dat = squeeze(res(vi,:,gi,:));
shadedErrorBar(Nvals, mean(dat), std(dat),'lineprops',{'-' 'color' c(gi,:)})
gi=2;
dat = squeeze(res(vi,:,gi,:));
shadedErrorBar(Nvals, mean(dat), std(dat),'lineprops',{'-' 'color' c(gi,:)})
gi=3;
dat = squeeze(res(vi,:,gi,:));
shadedErrorBar(Nvals, mean(dat), std(dat),'lineprops',{'-' 'color' c(gi,:)})

legend({'[0.5 0.5] \rho=0.2' ...
    '[0 0.5]\rho=0'...
    '[0.5 0.75] \rho=-0.2'},'location','northeast')
xlabel('Participants (N)')
ylabel('HPDI Width')
title('Within - 96% HPDI Width')
axis square
grid off

