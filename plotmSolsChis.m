function fig = plotmSolsChis(sols, styles, labels)
  fig = figure;
  set(fig, 'paperunits', 'centimeters', ...
    'papersize', [14 10], 'paperposition', [0 0 14 10]);

  for ind = 1:length(sols)
    sol = sols(ind);
    style = styles{ind};
    plot(sol.ts, sol.chis-1, style, 'linewidth', 1);
    hold on;
  end
  plot(sol.ts, sol.ts*0 + 0.25, 'k--', 'linewidth', 1);
  hold off;

  xlabel('t');
  ylabel('\chi(t)-1');
  ylim([0, max(sol.chis-1) * 1.1]);
  grid;

  legend(labels, 'location', 'southeast');
end

% print -color -depsc blasiusMultiFrs.eps
