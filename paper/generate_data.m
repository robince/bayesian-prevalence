function dat = generate_data(mu_g, sigma_b, sigma_w, Nsamp, Nsub)
% generate data from heirachical normal model
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
