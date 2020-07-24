gts = [0.25 0.25; 0.25 0.5; 0.25 0.75];
Nvals = 2:2:256;
% Nvals = 10;
a = 0.05;
b = 1;
Ngt = size(gts,1);
hpd = 0.96;

Nsamp = 1000;

parres = cell(1,length(Nvals));
parfor ni=1:length(Nvals)
    ni
    N = Nvals(ni);
    
    res = zeros(2,Nsamp,Ngt);
    % calcualte mean and s.d. for different ground truths
    for gi=1:Ngt
        theta1 = a + (b-a)*gts(gi,1);
        theta2 = a + (b-a)*gts(gi,2);
        
        k1 = binornd(N, theta1, [Nsamp 1]);
        k2 = binornd(N, theta2, [Nsamp 1]);
        
        for si=1:Nsamp
            [map, x, post, hpdi] = bayesprev_diff_between(k1(si),N,k2(si),N,hpd);
            res(1,si,gi) = map;
            res(2,si,gi) = hpdi(2) - hpdi(1);
        end
    end
    parres{ni} = res;
end

%%
res = cell2mat(reshape(parres,[1 1 1 length(Nvals)]));
save bayes_scale_between res Nvals gts Ngt hpd
