%% function fig = plotSol2(sol, bagSol, fign) - comparison between two solutions 
% (typically actual vs. Bagnoldian approximation)

function fig = plotSol2(sol, bagSol, fign)
  nContours = 60;

  if (nargin > 3 || nargin < 2)
      print_usage();
  end

  if (nargin == 3)
      fig = figure(fign);
  else
      fig = figure;
  end

  comparable = (sol.zg == bagSol.zg && sol.tg == bagSol.tg);

%% Contours of u(z,t)
%{
  subplot(2,2,1);
  contour(sol.tg, sol.zg, sol.ug, ...
    linspace(min (0, prctile(sol.ug(:), 0) ), max(sol.ug(:)), nContours), ...
    'LineWidth', 1.5 ...
  );

  hold on;
  contour(bagSol.tg, bagSol.zg, bagSol.ug, ...
    linspace(min (0, prctile(sol.ug(:), 5) ), max(sol.ug(:)), nContours), ...
    '--'     );
  hold off;

  view(2); colorbar('Location', 'SouthOutside'); colormap jet;
  title('u(z,t)');
  xlabel('t'); ylabel('z');
%}

%% Flow rate
  % subplot(2,2,2);
  subplot(2,1,1);
  dz = sol.zs(2) - sol.zs(1);
  plot(sol.ts, sol.qs, 'k-', ...
       bagSol.ts, bagSol.qs, 'k--');
  % grid; grid minor on;
  title('Flow rate q(t)');
  xlabel('t'); ylabel('q');

%% Shape factor
  % subplot(2,2,3);
  subplot(2,1,2);
  plot(sol.ts, sol.chis, 'k-', ...
       bagSol.ts, bagSol.chis, 'k--');
  % grid; grid minor on;
  title('Shape factor');
  xlabel('t'); ylabel('\chi');
  ylim([1.0, 1.35]);
  set(gca, 'ytick', 1.0:0.05:1.35);
  
%% Lyapunov function
%{
  if (comparable)
      subplot(2,2,4);
      zs = sol.zg(:, 1);
      nz = length(zs);
      ts = sol.tg(1, :);
      nt = length(ts);
      dzs = (zs(2:nz) - zs(1:nz-1)).'; dzs = [dzs, 0];
      Lyas = zeros(size(sol.ts));
      for ( tind = 1:length(sol.ts) )
          Lyas(tind) = (1/2)*integrate((sol.ug(:, tind) - bagSol.ug(:, tind)).^2, dz);
      end
      semilogy(ts, Lyas);
      % grid;
      title('Lyapunov');
      xlabel('t'); ylabel('L'); 
      % ylim([0, max(Lyas(ts > max(ts)/3))*sqrt(2)]);
      ylim([1e-6, inf]);
  end 
%}

end
