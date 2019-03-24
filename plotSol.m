%% function fig = plotSol(sol, fign)

function fig = plotSol(sol, alpha, fign)
  nContours = 60;

  if (nargin == 3)
      fig = figure(fign);
  else
      fig = figure;
  end

  if (nargin < 2)
      alpha = 0;
  end

  subplot(2,1,1);
  contour(sol.tg, sol.zg, sol.ug, ...
    linspace(min (0, prctile(sol.ug(:), 0) ), max(0, prctile(sol.ug(:), 100)), nContours) ...
  );
  view(2); colorbar('Location', 'SouthOutside'); colormap jet;

  title('u(z,t)');
  xlabel('t'); ylabel('z');

  subplot(2,1,2);
  plot(sol.ts, sol.chis-1, 'k-');
  title('\chi(t)-1');
  ylim([0, 0.4]);
  grid;

end
