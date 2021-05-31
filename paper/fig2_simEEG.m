% Ince, Paton, Kay and Schyns
% "Bayesian inference of population prevalence"
% biorxiv: https://doi.org/10.1101/2020.07.08.191106
%
% Figure 2: Simulated examples where Bayesian prevalence and second-level 
% t-tests diverge. EEG traces are simulated for 100 trials from 20 
% participants as white noise [N(0,1)] with an additive Gaussian activation 
% (σ = 20 ms) with amplitudes drawn from a uniform distribution. For each 
% simulation, mean traces are shown per participant (top left). A 
% second-level t-test is performed at each time point separately (blue 
% curve, bottom panel), dashed line shows the p=0.05 threshold, Bonferonni 
% corrected over time points. A within-participant t-test is performed at 
% each time point and for each participant separately (right hand panel); 
% the blue points show the maximum T-statistic over time points for each 
% participant, and the dashed line shows the p=0.05 Bonferonni corrected 
% threshold. Bottom right panel shows posterior population prevalence for 
% an effect in the analysis window. Black curves (bottom panel) show the 
% prevalence posterior at each time point (black line MAP, shaded area 
% 96% HPDI). A: An effect is simulated in all participants, with a peak 
% time uniformly distributed in the range 100-400 ms. B: An effect is 
% simulated in 10 participants, with a peak time uniformly distributed in 
% the range 200-275 ms.

Nsub = 20;
Ntrl = 100;

% Fs = 1000Hz (1ms bins)
Ntime = 600;

% 
% Panel A: all subjects have an effect but at different times
% 
% s = rng;
% save('figsubjectalingment','s')
load figsubjectalingment
rng(s);

% peak effects uniformly distributed between 100 and 500 ms
effect_range_width = 400;
effect_range_start = 100;
sub_effect_time = round(rand(Nsub,1)*effect_range_width + effect_range_start);
effect_size = 0.6;
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

stemlim = [3 14];
plot_results(dat, Nsub, Ntrl, Ntime, stemlim)


%
% Panel B: 50% of subjects have an effect
%
% s = rng;
% save('figsubjectprop','s')
load figsubjectprop
rng(s);

Nsubeff = 10;
Nsubnoeff = 10;
Nsub = Nsubeff + Nsubnoeff;

% peak effect uniformly distributed between 200 and 275 ms
effect_range_width = 75;
effect_range_start = 200;
sub_effect_time = round(rand(Nsub,1)*effect_range_width + effect_range_start);
effect_size = 0.6;
% generate data
dat = randn(Ntime,Ntrl,Nsub);
% add effect
effect_bins = 40;
effect_sd = 20;
for si=1:Nsubeff
    effect = normpdf(-effect_bins:effect_bins,0,effect_sd)';
    effect = effect./max(effect);
    % uniformly distributed random ampltiude on each trial 0-effect_size
    trial_amp = effect_size*rand(1,Ntrl);
    idx = sub_effect_time(si)-effect_bins:sub_effect_time(si)+effect_bins;
    dat(idx,:,si) = dat(idx,:,si) + effect.*trial_amp;
end

stemlim = [0 14];
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

% yyaxis right
% ylabel('Prev')
% ah = gca;
% ah.YAxis(2).Color = 'k';
% ah.YAxis(2).Label.Color = 'w';
% t = ah;

ah = subplot(3,3,[7 8]);
plot(group_t,'LineWidth',2)
cb = colorbar;
set(cb,'Vis','off')
ylim([-3 7])
xlabel('Time')
ylabel('T(19)')
yline(tinv(1-(0.05/Ntime),Nsub-1),'k--')
hold on
yyaxis right


x = 1:Ntime;
kt = sum(ind_sig,2);
y = zeros(1,Ntime);
lo = zeros(1,Ntime);
hi = zeros(1,Ntime);
for ti=1:Ntime
    y(ti) = bayesprev_map(kt(ti),Nsub);
    hp = bayesprev_hpdi(0.96,kt(ti),Nsub);
    lo(ti) = hp(1);
    hi(ti) = hp(2);
end
shadedErrorBar(x,y,[hi-y; y-lo],'lineprops',{'LineWidth',2,'Color','k'})
ylabel('Prevalence')

ah = gca;
ah.YAxis(1).Color = co(1,:);
ah.YAxis(2).Color = 'k';

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

