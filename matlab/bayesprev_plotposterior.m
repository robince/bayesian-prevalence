function h = bayesprev_plotposterior(k, n, a)
% Helper function to plot posterior distribution with MAP, 50% and 96% HPDI
% and 1st percentile of posterior.
%
% k : number of participants significant out of 
% n : total number of participants
% a : alpha value of within-participant test (default=0.05)


b = 1;
if nargin<3
    a = 0.05;
end

figure

% widths of HPDI indicators
oil = 3;
iil = 10;
% yaxis height of HPDI indicator
yp = 0.15;

x = linspace(0,1,100);
lw = 2;

co = get(gca,'ColorOrder');
cidx = get(gca,'ColorOrderIndex');

lh(1) = plot(x, bayesprev_posterior(x, k, n, a, b),'LineWidth',lw,'Color',co(cidx,:));
hold on
xmap = bayesprev_map(k,n, a, b);
pmap = bayesprev_posterior(xmap,k,n, a, b);
h96 = bayesprev_hpdi(0.96,k,n, a, b);

plot(xmap, yp,'.','MarkerSize',60,'Color',co(cidx,:));

plot([h96(1) h96(2)],[yp yp],'LineWidth',oil,'Color',co(cidx,:))
h50 = bayesprev_hpdi(0.5,k,n, a, b);
plot([h50(1) h50(2)],[yp yp],'LineWidth',iil,'Color',co(cidx,:))
% xline(xmap,'k')
box off

xlabel('Prevalence Proportion')
ylabel('Posterior Density')

title(sprintf('Posterior Prevalence from %d / %d at a=%.02f',k,n,a))

% 1st percetile
lb1 = bayesprev_bound(0.99,k,n,a);
xline(lb1,'Color',co(cidx,:))

fprintf(1,'\n %d / %d significant at a=%0.2f\n',k,n,a)
fprintf(1,'MAP [96%% HPDI]:  %.3f [%.3f %.3f]    [50%% HPDI]: [%.3f %.3f]\n',xmap,h96(1),h96(2),h50(1),h50(2))
fprintf(1,'1st percentile: %.3f\n\n', lb1)