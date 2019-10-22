% 
% Example of loading data from a CSV and applying Bayesian prevalence
% second level

% example_data.csv is a file with one binary value for each experimental
% unit (participant, neuron, voxel etc.) where 1 indicates the within-unit
% null hypothesis was rejected and 0 indicates it was not rejected

% load the data
sigdat = csvread('example_data.csv');

alpha = 0.05; % this specifies the alpha value used for the within-unit tests
Ntests = numel(sigdat); % number of tests (e.g. participants)
Nsigtests = sum(sigdat(:)); % number of significant tests


% example of Bayesian prevalence analyses on a single plot
figure
k = Nsigtests;
n = Ntests;
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
