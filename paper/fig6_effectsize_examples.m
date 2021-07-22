% Ince, Paton, Kay and Schyns
% "Bayesian inference of population prevalence"
% biorxiv: https://doi.org/10.1101/2020.07.08.191106
%
% Figure 6: Examples of different effect size prevalence curves with 
% similar p=0.05 prevalence. EEG traces are simulated for 100 trials from 
% 20 participants as white noise [N(0,1)] with an additive Gaussian 
% activation (σ = 20 ms) with amplitudes drawn from a uniform distribution. 
% For each simulation, mean traces are shown per participant (upper left 
% panel). A within-participant t-test is performed at each time point and 
% for each participant separately (right hand panel); the blue points show 
% the maximum T-statistic over time for each participant, and the dashed 
% line shows the p=0.05 Bonferroni corrected threshold. Lower right panel 
% shows posterior distribution of population prevalence for an effect in 
% the analysis window. Black curves (lower left panel) show the prevalence 
% (MAP, shaded area 96% HPDI) as a function of effect size threshold. A: 
% A weak early effect is simulated in all participants (peak time 
% uniformly distributed 100-150ms) B: In addition to the same early effect,
% a stronger, longer (σ = 40 ms) and more temporally variable later effect 
% is simulated in 10 participants, (peak times 250-450ms). C: Early events 
% are simulated with the same timing as A, but each participant has a 
% different maximum amplitude (participants ordered by effect size). All 
% three simulations have similar prevalence of p=0.05 effects, but show 
% differing patterns of prevalence over different effect size thresholds. 

% s = rng;
% save('figsimeeges','s')
load figsimeeges
rng(s);

Ntrl = 100;
% Fs = 1000Hz (1ms bins)
Ntime = 600;
Nsubeff = 10;
Nsubnoeff = 10;
Nsub = Nsubeff + Nsubnoeff;

% early weak effect in all participants
effect_range_width = 50;
effect_range_start = 100;
sub_effect_time = round(rand(Nsub,1)*effect_range_width + effect_range_start);
effect_size = 0.2;

% generate data
dat = randn(Ntime,Ntrl,Nsub);
% add effect
effect_bins = 40;
effect_sd = 20;
for si=1:Nsub
    effect = normpdf(-effect_bins:effect_bins,0,effect_sd)';
    effect = effect./max(effect);
    % uniformly distributed random ampltiude on each trial 0-effect_size
    trial_amp = effect_size*rand(1,Ntrl);
    idx = sub_effect_time(si)-effect_bins:sub_effect_time(si)+effect_bins;
    dat(idx,:,si) = dat(idx,:,si) + effect.*trial_amp;
end

stemlim = [0 20];
plot_results(dat, Nsub, Ntrl, Ntime, stemlim)

% add a later variable strong effect in some participants
effect_range_width = 200;
effect_range_start = 250;
sub_effect_time = round(rand(Nsub,1)*effect_range_width + effect_range_start);
% sub_effect_time = sort(sub_effect_time,'descend');
effect_size = 1;
% add effect
effect_bins = 80;
effect_sd = 40;
for si=1:Nsubeff
    effect = normpdf(-effect_bins:effect_bins,0,effect_sd)';
    effect = effect./max(effect);
    % uniformly distributed random ampltiude on each trial 0-effect_size
    trial_amp = effect_size*rand(1,Ntrl);
    idx = sub_effect_time(si)-effect_bins:sub_effect_time(si)+effect_bins;
    dat(idx,:,si) = dat(idx,:,si) + effect.*trial_amp;
end

plot_results(dat, Nsub, Ntrl, Ntime, stemlim)

% early weak effect in all participants, subject specific effect size
effect_range_width = 50;
effect_range_start = 100;
sub_effect_time = round(rand(Nsub,1)*effect_range_width + effect_range_start);
effect_size = linspace(1,0,Nsub);

% generate data
dat = randn(Ntime,Ntrl,Nsub);
% add effect
effect_bins = 40;
effect_sd = 20;
for si=1:Nsub
    effect = normpdf(-effect_bins:effect_bins,0,effect_sd)';
    effect = effect./max(effect);
    % uniformly distributed random ampltiude on each trial 0-effect_size
    trial_amp = effect_size(si)*rand(1,Ntrl);
    idx = sub_effect_time(si)-effect_bins:sub_effect_time(si)+effect_bins;
    dat(idx,:,si) = dat(idx,:,si) + effect.*trial_amp;
end
plot_results(dat, Nsub, Ntrl, Ntime, stemlim)



function plot_results(dat, Nsub, Ntrl, Ntime, stemlim)

% filter <30Hz as typical ERP analysis
Fs = 1000;
[b,a] = butter(3,30/(Fs/2),'low');
for si=1:Nsub
    for ti=1:Ntrl
        dat(:,ti,si) = filtfilt(b,a,dat(:,ti,si));
    end
end

% group results (t-test across participants at each time point)
group_t = zeros(Ntime,1);
group_p = zeros(Ntime,1);
for ti=1:Ntime
    [sig p ci stats] = ttest(squeeze(mean(dat(ti,:,:),2)));
    group_t(ti) = stats.tstat;
    group_p(ti) = p;
end

% group results (average over time points within each participant)
tavg = squeeze(mean(dat));
[sig p ci stats] = ttest(mean(tavg));
group_tavg_t = stats.tstat;
group_tavg_p = p;

% within participant results
ind_t = zeros(Ntime,Nsub);
ind_p = zeros(Ntime,Nsub);
for si=1:Nsub
    for ti=1:Ntime
        [sig p ci stats] = ttest(dat(ti,:,si));
        ind_t(ti,si) = stats.tstat;
        ind_p(ti,si) = p;
    end
end
ind_sig = ind_p<0.05/Ntime;

% % sig time points, unc and bonf
% [sum(group_p<0.05) sum(group_p<0.05/Ntime)]
% % sig subjects, unc and bonf 
% [sum(sum(ind_p<0.05)>0) sum(min(ind_p)<0.05/Ntime)]

% close all
figure
cm = flipud(cbrewer('div','RdBu',128));
co = get(gca,'ColorOrder');

imah = subplot(3,3,[1 2 4 5]);
imagesc(squeeze(mean(dat,2))')
colorbar
colormap(cm)
caxis([-1 1]*max(abs(caxis)))
set(gca,'YDir','normal')
yl = ylim;
% axis square
ylabel('Participant')

ah = subplot(3,3,[7 8]);
d.indt = max(ind_t);
d.Nsamp = size(dat,2);
d.Nsub = size(dat,3);
[es pmap hpdi] = prev_curve_onesided(d,1);
posbar = hpdi(2,:) - pmap;
negbar = pmap - hpdi(1,:);
shadedErrorBar(es,pmap,cat(1,posbar,negbar))
% xline(mu_g./(sigma_w./sqrt(Nsamp)),'b')
p = 0.05 ./ Ntime;
xline(tinv(1-p, d.Nsamp-1),'r');
% xline(tinv(p/2,d.Nsamp-1),'r');
ylim([0 1])
xlim([0 20])
xlabel('Threshold T(99)')
ylabel('Prevalence (> Threshold)')
cb = colorbar;
set(cb,'Vis','off')

subplot(3,3,[3 6])
% plot(max(ind_t,[],1),'s','MarkerSize',10,'LineWidth',2)
stem(max(ind_t,[],1),'filled','LineWidth',1);
set(gca,'view',[90 -90])
ylim(stemlim)
xlim(yl)
ylabel('T(99)')
h = hline(tinv(1-(0.05/Ntime),Ntrl-1),'k--');
% set(h,'LineWidth',1.5)

subplot(3,3,9)
k = sum((min(ind_p)<0.05/Ntime));
i=1;
oil = 5;
iil = 15;
hy = 0.3;
a = 0.05;
b = 1;
co = get(gca,'ColorOrder');
x = linspace(0,1,100);
lw = 2;

lh(1) = plot(x, bayesprev_posterior(x, k, Nsub, a, b),'Color','k','LineWidth',lw);
hold on

xmap = bayesprev_map(k, Nsub, a, b);
pmap = bayesprev_posterior(xmap, k, Nsub, a, b);    
h = bayesprev_hpdi(0.96,k, Nsub, a, b);

% yp = pmap;
yp = 0.5;
yp = 0.25;
c = [0 0 0 0.4];
plot(xmap, yp,'.','MarkerSize',20,'Color','k');
plot([h(1) h(2)],[yp yp],'Color',c,'LineWidth',oil)
h = bayesprev_hpdi(0.5,k,Nsub, a, b);
plot([h(1) h(2)],[yp yp],'Color',c,'LineWidth',iil)
% xline(xmap,'k')
box off
xlabel('Population Prevalence')
ylabel('Posterior Density')

end



