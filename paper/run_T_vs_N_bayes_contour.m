
sigma_w = 10;
sigma_b = 2;
mu_g = 1;
% Nvals = 2:2:200;
% kvals = 2:2:500;

% Nvals = 2.^[1:8];
% kvals = 2.^[1:9];

Nvals = [2 4 8 16 32 64 100 150 200 250];
kvals = [2 4 8 16 32 64 128 200 250 300 350 400 450 500];

NN = length(Nvals);
Nk = length(kvals);
Nperm = 1000;

gmap = zeros(Nperm,Nk,NN);
glb = zeros(Nperm,Nk,NN);

tic
parfor ni=1:NN
    ni
    for ki=1:Nk
        Nsub = Nvals(ni);
        Nsamp = kvals(ki);
        for pi=1:Nperm
            submeanstrue = normrnd(mu_g, sigma_b, [Nsub 1]);
            indsig = zeros(1,Nsub);
            for si=1:Nsub
                dat = normrnd(submeanstrue(si), sigma_w, [Nsamp 1]);
                indsig(si) = ttest(dat);
            end
            k = sum(indsig);
            gmap(pi,ki,ni) = bayesprev_map(k,Nsub);
            glb(pi,ki,ni) = bayesprev_bound(0.95,k,Nsub);
        end
    end
end
toc


save prevbayes_normal Nvals kvals gmap glb