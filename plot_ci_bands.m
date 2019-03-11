function h = plot_ci_bands(x,Y, color, plotfun)

%yb = std(Y,1,2);
yb_upper = prctile(Y,75,2);
yb_lower = prctile(Y,25,2);
yavg = median(Y,2);


h = plotfun(x, yavg, 'color', color, 'LineWidth',2);
hold on
fill_between(x, yb_upper, yb_lower, color, 'LineStyle','none')
alpha(.1)
%%
