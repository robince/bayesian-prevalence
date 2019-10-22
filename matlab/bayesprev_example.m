% Example of how to use Bayesian prevalence functions
%
% 1.   Simulate or load within-participant raw experimental data
% 2.   LEVEL 1: Apply statstical test at the individual level
% 3.   LEVEL 2: Apply Bayesian Prevalence to the outcomes of Level 1

%
% 1.   Simulate or load within-participant raw experimental data
%

% 1.1. Simulate within-participant raw experimental data
Nsub = 20; % number of particpants
Nsamp = 100; % trials/samples per participant
sigma_w = 10; % within-participant SD
sigma_b = 2; % between-participant SD
mu_g = 1; % population mean

% per participant mean drawn from population normal distribution
submeanstrue = normrnd(mu_g, sigma_b, [Nsub 1]);
% rawdat holds trial data for each participnt 
rawdat = zeros(Nsamp, Nsub);
for si=1:Nsub
    % generate trials for each participant
    rawdat(:,si) = normrnd(submeanstrue(si), sigma_w, [Nsamp 1]);
end

% 1.2.Load within-participant raw experimental data
% Load your own data into the variable rawdat with dimensions [Nsamp Nsub],
% setting Nsamp and Nsub accordingly.

%
% 2.  LEVEL 1
%

% 2.1. Within-participant statistical test
% The loop performs within-participant statistical test. Here, a t-test for
% non-zero mean which is the simplest statistical test. In general, any
% statistical test can be used at Level 1.

indsig = zeros(1,Nsub);
for si=1:Nsub
    % within-participant t-test significance
    [indsig(si) p ci stats] = ttest(rawdat(:,si));
end
% the binary variable indsig indicates whether the within-participant
% t-test is, or not, significant for each participant (1 entry for each
% participant)

% 2.2. Loading within-participant statistical test.
% You can also load your own within-participant statistical test results here.
% Load binary results into binary indsig vector with one entry per participant.
% See also example_csv.m for an example.

%
% 3.  LEVEL 2
%

% Bayesian prevalence inference is performed with three numbers: 
% k, the number of significant participants (e.g. sum of binary indicator
% variable)
% n, the number of participants in the sample
% alpha, the false positive rate
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
% add 50% HPDI
h = bayesprev_hpdi(0.5,k,n,alpha);
plot([h(1) h(2)],[pmap pmap],'Color',co(ci,:),'LineWidth',iil)

xlabel('Population prevalence proportion')
ylabel('Posterior density')
