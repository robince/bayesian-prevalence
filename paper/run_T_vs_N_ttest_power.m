sigma_w = 20;
sigma_b = 2;
mu_g = 1;

Nvals = 2:2:200;
kvals = 2:2:500;

NN = length(Nvals);
Nk = length(kvals);

tpow = zeros(Nk,NN);
parfor ni=1:NN
    for ki=1:Nk
        sigma_g = sqrt(sigma_b.^2 + ((sigma_w).^2)/kvals(ki));
        tpow(ki,ni) = sampsizepwr('t',[0 sigma_g],mu_g,[],Nvals(ni));
    end
end

%%
save tpow tpow Nvals kvals
