%% function fig = plotSol(sol, fign)

function fig = plotSolChi(sol, alpha, fign)
  if (nargin == 3)
      fig = figure(fign);
  else
      fig = figure;
  end

  if (nargin < 2)
      alpha = 0;
  end

  loglog(sol.ts, (sol.chis-1).^(4/3), 'k-', ...
       sol.ts, sol.ts*0 + .25^(4/3), 'k--');
  xlabel('t');
  ylabel('\chi(t)-1');
  ylim([1e-3, max(sol.chis-1) * 1.1]);
  grid;
end

