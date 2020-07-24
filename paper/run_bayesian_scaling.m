
Nvals = 2:2:256;
a = 0.05;
b = 1;

gts = [0.25 0.5 0.75];
Ngt = length(gts);
hpd = 0.96;

parres = cell(1,length(Nvals));
    
parfor ni=1:length(Nvals)
    ni
    N = Nvals(ni);
    k = 0:N;
    hpdiwidthk = zeros(1,N+1);
    lboundk = zeros(1,N+1);
    for ki=1:N+1
        hpdi = bayesprev_hpdi(hpd,k(ki),N);
        hpdiwidthk(ki) = hpdi(2)-hpdi(1);
        lboundk(ki) = bayesprev_bound(0.95, k(ki), N);
    end
    
    res = zeros(2,2,Ngt);
    % calcualte mean and s.d. for different ground truths
    for gi=1:Ngt
        theta = a + (b-a)*gts(gi);
        pk = binopdf(k, N, theta);
        mu = sum(pk.*hpdiwidthk);
        res(1,1,gi) = mu;
        res(2,1,gi) = sqrt(sum(pk.*(hpdiwidthk-mu).^2));
        mu = sum(pk.*lboundk);
        res(1,2,gi) = mu;
        res(2,2,gi) = sqrt(sum(pk.*(lboundk-mu).^2));
    end
    parres{ni} = res;
end

%%
res = cell2mat(reshape(parres,[1 1 1 length(Nvals)]));
save bayes_scale res Nvals gts Ngt hpd
