%% function fig = plotSol(sol, fign)

function ax = plotSolChi(sol, style, ax)
  if (nargin < 2 || strcmp(style, ''))
      style = 'k';
  end
  if (nargin < 3)
      ax = gca;
  end

  plot(ax, sol.ts, sol.chis-1, style, 'linewidth', 2, ...
       sol.ts, sol.ts*0 + 0.25, 'k--', 'linewidth', 2);
  xlabel('t');
  ylabel('\chi(t)-1');
  ylim([0, max(sol.chis-1) * 1.1]);
  grid;
  
  ax = gca;
end

