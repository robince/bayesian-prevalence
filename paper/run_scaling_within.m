
gts = [0.5 0.5 0.2; 0 0.5 0; 0.5 0.75 -0.2];

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
        g1 = gts(gi,1);
        g2 = gts(gi,2);
        g11 = g1*g2 + gts(gi,3)*sqrt(g1*(1-g1)*g2*(1-g2));
        g10 = g1 - g11;
        g01 = g2 - g11;
        g00 = 1 - g11 - g10 - g01;
        
        the11 = (b^2)*g11 + a*b*g10 + a*b*g01 + a*a*g00;
        the10 = a + (b-a)*g1 - the11;
        the01 = a + (b-a)*g2 - the11;
        the00 = 1 - the11 - the10 - the01;
        theta = [the11 the10 the01 the00];
        
        for si=1:Nsamp
            k = mnrnd(N, theta);
            [map, x, post, hpdi] = bayesprev_diff_within(k(1),k(2),k(3),N,hpd);
            res(1,si,gi) = map;
            res(2,si,gi) = hpdi(2) - hpdi(1);
        end
    end
    parres{ni} = res;
end

%%
res = cell2mat(reshape(parres,[1 1 1 length(Nvals)]));
save bayes_scale_within res Nvals gts Ngt hpd
