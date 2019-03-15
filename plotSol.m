%% function fig = plotSol(sol, bagSol, fign)

function fig = plotSol(sol, bagSol, fign)
  nContours = 60;

  if (nargin > 3)
      print_usage();
  end

  if (nargin == 3)
      fig = figure(fign);
  else
      fig = figure;
  end

  if (nargin == 1)
    subplot(2,1,1);
    contour(sol.tg, sol.zg, sol.ug, ...
      linspace(min (0, prctile(sol.ug(:), 5) ), max(0, prctile(sol.ug(:), 95)), nContours) ...
    );
    view(2); colorbar('Location', 'SouthOutside'); colormap jet;

    % mesh(sol.tg, sol.zg, real(sol.ug));
    title('u(z,t)');
    xlabel('t'); ylabel('z');

    subplot(2,1,2);
    dz = sol.zs(2) - sol.zs(1);
    chis = 1/dz * sum( sol.ug.^2 , 1) ./ ( sum(sol.ug, 1).^2 );
    plot(sol.ts, chis, 'k-');
    ylim([1, 1.6]);
  end

%% Contours of u(z,t)
  if (nargin >= 2)
    comparable = (sol.zg == bagSol.zg && sol.tg == bagSol.tg);
    if (comparable)
        subplot(2,2,1);
    else
        subplot(3,1,1);
    end
    contour(sol.tg, sol.zg, sol.ug, ...
      linspace(min (0, prctile(sol.ug(:), 5) ), max(sol.ug(:)), nContours), ...
      'LineWidth', 1.5 ...
    );
    %{
    hold on;
    contour(bagSol.tg, bagSol.zg, bagSol.ug, ...
      linspace(min (0, prctile(sol.ug(:), 5) ), max(sol.ug(:)), nContours), ...
      ':'     );
    hold off;
    %}
    % mesh(sol.tg, sol.zg, real(sol.ug));
    view(2); colorbar('Location', 'SouthOutside'); colormap jet;
    title('u(z,t)');
    xlabel('t'); ylabel('z');

    %{
%% Difference-squared
    if (comparable)
        subplot(2,2,2);
    else
        subplot(3,1,2);
    end
    diffsq = (sol.ug - bagSol.ug).^2;
    contour(sol.tg, sol.zg, diffsq, ...
      linspace(prctile(diffsq(:), 5), prctile(diffsq(:), 95), nContours) );
    view(2); colorbar('Location', 'SouthOutside'); colormap jet;
    title('(u - B)^2');
    xlabel('t'); ylabel('z');
    %}

%% Flow rate
    subplot(2,2,2);
    dz = sol.zs(2) - sol.zs(1);
    qs = zeros(size(sol.ts));
    for ( tind = 1:length(sol.ts) )
        qs(tind) = integrate(sol.ug(:, tind), dz);
    end
    plot(sol.ts, qs, 'k-');
    grid;
    grid minor on;
    title('Flow rate q(t)');
    xlabel('t'); ylabel('q');



%% Shape factor
    subplot(2,2,3);
    dz = sol.zs(2) - sol.zs(1);
    chis = zeros(size(sol.ts));
    for ( tind = 1:length(sol.ts) )
        chis(tind) = integrate(sol.ug(:, tind).^2, dz)  ...
                        / integrate(sol.ug(:, tind), dz)^2; 
    end
    plot(sol.ts, chis, 'k-');
    grid;
    grid minor on;
    title('Shape factor');
    xlabel('t'); ylabel('\chi');
    ylim([1.0, 1.35]);
    set(gca, 'ytick', 1.0:0.05:1.35);
    
%% Lyapunov function
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
        plot(ts, Lyas);
        title('Lyapunov');
        xlabel('t'); ylabel('L'); 
        % ylim([0, max(Lyas(ts > max(ts)/3))*sqrt(2)]);
        ylim([0, inf]);
    end 

  end
end
