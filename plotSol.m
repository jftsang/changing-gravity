%% function fig = plotSol(sol, fign)

function fig = plotSol(sol, alpha, fign)

  if (nargin == 3)
      fig = figure(fign);
  else
      fig = figure;
  end

  if (nargin < 2)
      alpha = 0;
  end

  nContours = 60;
  contourVals = linspace(min(0, prctile(sol.ug(:), 0) ), ...
                         max(0, prctile(sol.ug(:), 100)), nContours);
  subplot(2,1,1);
  contour(sol.tg(:,2:end), sol.zg(:,2:end), ... 
    sol.ug(:, 2:end), ...
    contourVals ...
  );
  view(2); colorbar('Location', 'SouthOutside'); colormap jet;
  title('u(z,t)');
  xlabel('t'); ylabel('z');
  ylim([0, .15]);

  subplot(2,1,2);
  plot(sol.ts, (sol.chis-1), 'k-', ...
       sol.ts, sol.ts*0 + .25, 'k--');
  xlabel('t');
  ylabel('\chi(t)-1');
  ylim([0, max(sol.chis-1) * 1.1]);
  grid;
end
