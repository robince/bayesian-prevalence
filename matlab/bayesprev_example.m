% Example of how to use Bayesian prevalence functions

% simulate some data

Nsub = 20; % number of particpants
Nsamp = 100; % trials/samples per participant
sigma_w = 10; % within-participant SD
sigma_b = 2; % between-participant SD

% generate within-participant data
submeanstrue = normrnd(mu_g, sigma_b, [Nsub 1]);
rawdat = zeros(Nsamp, Nsub);
for si=1:Nsub
    rawdat(:,si) = normrnd(submeanstrue(si), sigma_w, [Nsamp 1]);
end

% perform within-participant statistical test (here t-test)
indsig = zeros(1,Nsub);
for si=1:Nsub
    % within-participant t-test significance
    [indsig(si) p ci stats] = ttest(rawdat(:,si));
end

% within-participant test results for prevalence inference 
k = sum(indsig);
n = Nsub;
alpha = 0.05; % default value see 'help ttest'

% plot posterior distribution of population prevalence
figure
co = get(gca,'ColorOrder'); ci=1;
hold on

x = linspace(0,1,100);
posterior = bayesprev_posterior(x,k,n,alpha);
plot(x, posterior,'Color',co(ci,:));

% add MAP as a point
xmap = bayesprev_map(k,n,alpha);
pmap = bayesprev_posterior(xmap,k,n,alpha);
plot(xmap, pmap,'.','MarkerSize',20,'Color',co(ci,:));

% add lower bound as a vertical line
bound = bayesprev_bound(0.95,k,n,alpha);
line([bound bound], [0 bayesprev_posterior(bound,k,n,alpha)],'Color',co(ci,:),'LineStyle',':')

% add 95% HPDI
oil = 2;
iil = 4;
h = bayesprev_hpdi(0.95,k,n,alpha);
plot([h(1) h(2)],[pmap pmap],'Color',co(ci,:),'LineWidth',oil)
h = bayesprev_hpdi(0.5,k,Nsub, a, b);
plot([h(1) h(2)],[pmap pmap],'Color',co(ci,:),'LineWidth',iil)

xlabel('Population prevalence proportion')
ylabel('Posterior density')
