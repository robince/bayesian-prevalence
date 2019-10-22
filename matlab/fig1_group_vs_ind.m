%
% Figure 1: Population vs individual inference. For each simulation we sample ?
% = 50 individual participant mean effects from a normal distribution with
% population mean ? (A,B: ? = 0; C,D: ? = 1) and between-participant
% standard deviation ?_b = 2. Within each participant, ? trials (A,C: ? = 20;
% B,D: ? = 500) are drawn from a normal distribution with the
% participant-specific mean and a common within-participant standard deviation 
% ?_w = 10.Orange and blue indicate, respectively, exceeding or not exceeding a
% ?=0.05 threshold for a t-test at the population level (on the within-participant
% means, population normal density curves) or at the individual participant
% level (individual sample means +/- s.e.m.). E: Bayesian posterior population
% prevalence distributions for the 4 simulated data sets. Points show Bayesian
% maximum a posteriori estimate. Thick and thin horizontal lines indicate 50%
% and 96% highest posterior density intervals respectively.

% Simulations inspired by 
% Baker, Vilidaite, Lygo, Smith, Flack, Gouws and Andrews
% "Power contours: optimising sample size and precision in experimental
% psychology and human neuroscience"

x = [];

figure

Nsub = 50;
sigma_w = 10;
sigma_b = 2;

% s = rng;
% save('currfig1seed','s')
load currfig1seed
rng(s);

% Generate data from heirachical normal distribution
% and plot per-participant means with standard deviations 
% together with population distribution
subplot(3,2,1)
Nsamp = 20;
mu_g = 0;
sigma_g = sqrt(sigma_b.^2 + ((sigma_w).^2)/Nsamp);
datA = generate_data(mu_g, sigma_b, sigma_w, Nsamp, Nsub);
datA.sigma_g = sigma_g;
plot_data(mu_g, sigma_g, datA, datA.groupsig+1);

subplot(3,2,2)
Nsamp = 500;
mu_g = 0;
sigma_g = sqrt(sigma_b.^2 + ((sigma_w).^2)/Nsamp);
datB = generate_data(mu_g, sigma_b, sigma_w, Nsamp, Nsub);
datB.sigma_g = sigma_g;
plot_data(mu_g, sigma_g, datB, datB.groupsig+1);

subplot(3,2,3)
Nsamp = 20;
mu_g = 1;
sigma_g = sqrt(sigma_b.^2 + ((sigma_w).^2)/Nsamp);
datC = generate_data(mu_g, sigma_b, sigma_w, Nsamp, Nsub);
datC.sigma_g = sigma_g;
plot_data(mu_g, sigma_g, datC, datC.groupsig+1);

subplot(3,2,4)
Nsamp = 500;
mu_g = 1;
sigma_g = sqrt(sigma_b.^2 + ((sigma_w).^2)/Nsamp);
datD = generate_data(mu_g, sigma_b, sigma_w, Nsamp, Nsub);
datD.sigma_g = sigma_g;
plot_data(mu_g, sigma_g, datD, datD.groupsig+1);



subplot(3,1,3);
x = linspace(0,1,200);
co = get(gca,'ColorOrder');

oil = 2;
iil = 4;
a = 0.05;
lh = [];

k = sum(datA.indsig);i=3;hy = 0.3;
b = 1;
dat = datA;
lh(1) = plot(x, bayesprev_posterior(x, k, Nsub, a, b),'Color',co(i,:));
hold on
xmap = bayesprev_map(k,Nsub, a, b);
pmap = bayesprev_posterior(xmap,k,Nsub, a, b);
plot(xmap, pmap,'.','MarkerSize',20,'Color',co(i,:));
h = bayesprev_hpdi(0.96,k,Nsub, a, b);
plot([h(1) h(2)],[pmap pmap],'Color',co(i,:),'LineWidth',oil)
h = bayesprev_hpdi(0.5,k,Nsub, a, b);
plot([h(1) h(2)],[pmap pmap],'Color',co(i,:),'LineWidth',iil)

k = sum(datB.indsig);i=4;hy = 0.5;
b = 1;
dat = datB;
lh(2) = plot(x, bayesprev_posterior(x, k, Nsub, a, b),'Color',co(i,:));
hold on
xmap = bayesprev_map(k,Nsub, a, b);
pmap = bayesprev_posterior(xmap,k,Nsub, a, b);
plot(xmap, pmap,'.','MarkerSize',20,'Color',co(i,:));
h = bayesprev_hpdi(0.96,k,Nsub, a, b);
plot([h(1) h(2)],[pmap pmap],'Color',co(i,:),'LineWidth',oil)
h = bayesprev_hpdi(0.5,k,Nsub, a, b);
plot([h(1) h(2)],[pmap pmap],'Color',co(i,:),'LineWidth',iil)

k = sum(datC.indsig);i=5;hy = 0.5;
b = 1;
dat = datC;
lh(3) = plot(x, bayesprev_posterior(x, k, Nsub, a, b),'Color',co(i,:));
hold on
xmap = bayesprev_map(k,Nsub, a, b);
pmap = bayesprev_posterior(xmap,k,Nsub, a, b);
plot(xmap, pmap,'.','MarkerSize',20,'Color',co(i,:));
h = bayesprev_hpdi(0.96,k,Nsub, a, b);
plot([h(1) h(2)],[pmap pmap],'Color',co(i,:),'LineWidth',oil)
h = bayesprev_hpdi(0.5,k,Nsub, a, b);
plot([h(1) h(2)],[pmap pmap],'Color',co(i,:),'LineWidth',iil)

k = sum(datD.indsig);i=6;hy = 0.3;
b = 1;
dat = datD;
lh(4) = plot(x, bayesprev_posterior(x, k, Nsub, a, b),'Color',co(i,:));
hold on
xmap = bayesprev_map(k,Nsub, a, b);
pmap = bayesprev_posterior(xmap,k,Nsub, a, b);
plot(xmap, pmap,'.','MarkerSize',20,'Color',co(i,:));
h = bayesprev_hpdi(0.96,k,Nsub, a, b);
plot([h(1) h(2)],[pmap pmap],'Color',co(i,:),'LineWidth',oil)
h = bayesprev_hpdi(0.5,k,Nsub, a, b);
plot([h(1) h(2)],[pmap pmap],'Color',co(i,:),'LineWidth',iil)

xlim([0 1])

legend(lh,{'A','B','C','D'})

% generate data from heirachical normal model
function dat = generate_data(mu_g, sigma_b, sigma_w, Nsamp, Nsub)
% mu_g - ground truth group mean
% sigma_b - between participant standard deviation
% sigma_w - within participant standard deviation
% Nsamp - number of trials per participant
% Nsub - number of participants

% generate individual subject means from population normal distribution
submeanstrue = normrnd(mu_g, sigma_b, [Nsub 1]);
rawdat = zeros(Nsamp, Nsub);
dat.indsig = false(1,Nsub);
dat.indt = zeros(1,Nsub);
for si=1:Nsub
    % generate within-participant data
    rawdat(:,si) = normrnd(submeanstrue(si), sigma_w, [Nsamp 1]);
    % within-participant t-test significance
    [dat.indsig(si) p ci stats] = ttest(rawdat(:,si));
    % within-participant t-score
    dat.indt(si) = stats.tstat;
end
% within-participant mean
dat.submeans = mean(rawdat,1);
dat.subsem = std(rawdat,[],1) ./ sqrt(Nsamp);
% second level t-test on within-participant means
[h p ci stats] = ttest(mean(rawdat,1));
% population level t-test significance
dat.groupsig = h;
dat.groupp = p;
dat.Nsub = Nsub;
dat.Nsamp = Nsamp;
dat.tdf = stats.df;
% population level t-score
dat.t = stats.tstat;
end

% plot group and individual data
function plot_data(mu_g, sigma_g, d, gc)

x = -15:0.1:15;
y = normpdf(x,mu_g,sigma_g);
co = get(gca,'ColorOrder');
% plot population model (normal) distribution
plot(x,y,'Color',co(gc,:))
hold on
ypos = unifrnd(zeros(1,d.Nsub)+0.0005, 0.95*normpdf(d.submeans, mu_g, sigma_g));
% plot(submeans, ypos, '.' , 'Color',co(1,:))
% plot significant within-participant means in orange
errorbar(d.submeans(d.indsig),ypos(d.indsig),d.subsem(d.indsig),'horizontal','.',...
                                      'capsize',0,...
                                      'Color',co(2,:),...
                                      'MarkerSize',10)
% plot non-significant within-participant means in blue 
errorbar(d.submeans(~d.indsig),ypos(~d.indsig),d.subsem(~d.indsig),'horizontal','.',...
                                      'capsize',0,...
                                      'Color',co(1,:),...
                                      'MarkerSize',10)
xline(0,'k:');
if mu_g~=0
    xline(mu_g,'k--');
end
axis square
xlim([-15 15])
title(sprintf('%d trials, %d/50 sig, group t(%d)=%.2f p=%.3f',d.Nsamp,sum(d.indsig),d.tdf,d.t,d.groupp))

end
