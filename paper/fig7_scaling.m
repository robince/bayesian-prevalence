% Ince, Paton, Kay and Schyns
% "Bayesian inference of population prevalence"
% biorxiv: https://doi.org/10.1101/2020.07.08.191106
%
% Figure 7: Characterisation of Bayesian prevalence inference.  
% A,B,C: We consider the binomial model of within-participant testing for 
% three ground truth population proportions: 25%, 50% and 75% (blue, 
% orange, yellow, respectively). We show how A: the Bayesian MAP estimate, 
% B: 95% Bayesian lower bound and C: 96% HPDI width, scale with the number 
% of participants. Lines show theoretical expectation, coloured regions 
% show +/- 1 s.d.. D,E,F: We consider the population model from 
% Figure 1C,D (μ=1). D: Power contours for the population inference using 
% a t-test (Baker et al. 2020). Colour scale shows statistical power 
% (probability of rejecting the null hypothesis). E: Contours of average 
% Bayesian MAP estimate for γ. Colour scale shows MAP prevalence 
% proportion.  F: Contours of average 95% Bayesian lower bound for γ. 
% Colour scale shows lower bound prevalence. From the prevalence 
% perspective, the number of trials obtained per participant has a larger 
% effect on the resulting population inference than does the number of 
% participants.


figure

%
% Panel D, T-test power scaling with T and N
%
% tpow.mat from run_T_vs_N_ttest_power.m
subplot(2,3,4)
load tpow

cmap = cbrewer('seq','Oranges',100);
set(gca,'YDir','normal')
cb = colorbar;
colormap(cmap);
caxis([-0.5 1])
set(cb,'YLim',[0 1])
hold on 
st = 'on';
contour(Nvals,kvals,tpow,[0.1:0.1:0.9 0.99 0.9999],'ShowText',st)
xlabel('Participants (N)')
ylabel('Trials per participant (k)')
title('t-test power')
axis square

%
% Panel E: Bayesian prevalence MAP
%
% prevbayes_normal.mat from run_T_vs_N_bayes_contour.m
subplot(2,3,5)
load prevbayes_normal

cmap = cbrewer('seq','Oranges',100);
% cmap = flipud(cmap);
% cmap = cmap(300:400,:);
% imagesc(kvals,Nvals,tpow)
set(gca,'YDir','normal')
cb = colorbar;
colormap(cmap);
caxis([-0.5 1])
set(cb,'YLim',[0 1])
hold on 
st = 'on';
t = squeeze(mean(gmap,1));
xGrid = repmat(Nvals,numel(kvals),1);
yGrid = repmat(kvals',1,numel(Nvals));
[xQuery, yQuery] = meshgrid(1:250,1:500);
vq = interp2(xGrid,yGrid,t,xQuery,yQuery,'makima');
contour(xQuery,yQuery,vq,[0.1:0.1:0.6],'ShowText',st)
xlabel('Participants (N)')
ylabel('Trials per participant (k)')
title('MAP g')
axis square

%
% Panel F: Bayesian prevalence lower bound
%
% prevbayes_normal.mat from run_T_vs_N_bayes_contour.m
subplot(2,3,6)

cmap = cbrewer('seq','Oranges',100);
set(gca,'YDir','normal')
cb = colorbar;
colormap(cmap);
caxis([-0.5 1])
set(cb,'YLim',[0 1])
hold on 
st = 'on';
t = squeeze(mean(glb,1));
xGrid = repmat(Nvals,numel(kvals),1);
yGrid = repmat(kvals',1,numel(Nvals));
[xQuery, yQuery] = meshgrid(1:250,1:500);
vq = interp2(xGrid,yGrid,t,xQuery,yQuery,'makima');
contour(xQuery,yQuery,vq,[0.1:0.1:0.6],'ShowText',st)
xlabel('Participants (N)')
ylabel('Trials per participant (k)')
title('95% lower bound')
axis square


%
% Panel A: scaling of MAP with N
%
ax2 = [];
ax2(1) = subplot(2,3,1);
hold all

Nvals = 2:2:256;
a = 0.05;
b = 1;
c = get(0, 'DefaultAxesColorOrder');

gt = 0.25;
theta = a + (b-a)*gt;
dat = zeros(length(Nvals),2);
for ni=1:length(Nvals)
    N = Nvals(ni);
    k = 0:N;
    vk = zeros(1,N+1);
    for ki=1:N+1
        vk(ki) = bayesprev_map(k(ki),N);
    end
    pk = binopdf(k, N, theta);
    mu = sum(pk.*vk); % mean
    sigma = sqrt(sum(pk.*(vk-mu).^2));
    dat(ni,1) = mu;
    dat(ni,2) = sigma;
end
shadedErrorBar(Nvals, dat(:,1), dat(:,2),'lineprops',{'-' 'color' c(1,:)})

gt = 0.5;
theta = a + (b-a)*gt;
dat = zeros(length(Nvals),2);
for ni=1:length(Nvals)
    N = Nvals(ni);
    k = 0:N;
    vk = zeros(1,N+1);
    for ki=1:N+1
        vk(ki) = bayesprev_map(k(ki),N);
    end
    pk = binopdf(k, N, theta);
    mu = sum(pk.*vk); % mean
    sigma = sqrt(sum(pk.*(vk-mu).^2));
    dat(ni,1) = mu;
    dat(ni,2) = sigma;
end
shadedErrorBar(Nvals, dat(:,1), dat(:,2),'lineprops',{'-' 'color' c(2,:)})

gt = 0.75;
theta = a + (b-a)*gt;
dat = zeros(length(Nvals),2);
for ni=1:length(Nvals)
    N = Nvals(ni);
    k = 0:N;
    vk = zeros(1,N+1);
    for ki=1:N+1
        vk(ki) = bayesprev_map(k(ki),N);
    end
    pk = binopdf(k, N, theta);
    mu = sum(pk.*vk); % mean
    sigma = sqrt(sum(pk.*(vk-mu).^2));
    dat(ni,1) = mu;
    dat(ni,2) = sigma;
end
shadedErrorBar(Nvals, dat(:,1), dat(:,2),'lineprops',{'-' 'color' c(3,:)})

% legend({'25%' '50%' '75%'},'location','southeast')
xlabel('Participants (N)')
ylabel('\gamma')
title('Bayesian MAP')
axis square

%
% Panel B: scaling of 95% lower bound with N
%
% bayes_scale.mat from run_bayesian_scaling.m
ax2(2) = subplot(2,3,2);
hold all
load bayes_scale

plt = 2; 
gi=1;
shadedErrorBar(Nvals, res(1,plt,gi,:), res(2,plt,gi,:),'lineprops',{'-' 'color' c(gi,:)})
gi=2;
shadedErrorBar(Nvals, res(1,plt,gi,:), res(2,plt,gi,:),'lineprops',{'-' 'color' c(gi,:)})
gi=3;
shadedErrorBar(Nvals, res(1,plt,gi,:), res(2,plt,gi,:),'lineprops',{'-' 'color' c(gi,:)})

% legend({'25%' '50%' '75%'},'location','southeast')
xlabel('Participants (N)')
ylabel('95% lower bound')
title('Bayesian 95% lower bound')
axis square


%
% Panel C: scaling of 96% HPDI width with N
%
% bayes_scale.mat from run_bayesian_scaling.m
ax2(3) = subplot(2,3,3);
hold all
load bayes_scale

plt = 1; 
gi=1;
shadedErrorBar(Nvals, res(1,plt,gi,:), res(2,plt,gi,:),'lineprops',{'-' 'color' c(gi,:)})
gi=2;
shadedErrorBar(Nvals, res(1,plt,gi,:), res(2,plt,gi,:),'lineprops',{'-' 'color' c(gi,:)})
gi=3;
shadedErrorBar(Nvals, res(1,plt,gi,:), res(2,plt,gi,:),'lineprops',{'-' 'color' c(gi,:)})

legend({'25%' '50%' '75%'},'location','northeast')
xlabel('Participants (N)')
ylabel('HPDI width')
title('96% HPDI width')
axis square

set(ax2,'xlim',[0 max(Nvals)])
set(ax2,'ylim',[-0.05 1.1])
