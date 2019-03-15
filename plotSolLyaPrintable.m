%% function fig = plotSolLyaPrintable(sol, bagSol)

function fig = plotSolLyaPrintable(sol, bagSol)
  nContours = 60;

  fig = figure( ...
            'paperunits', 'centimeters', ...
            'papersize', [18 6], ...
            'paperposition', [1 0.5 16 5], ...
            'outerposition', [0 0 0 0] ...
        );

    dz = sol.zs(2) - sol.zs(1);
    Lyas = zeros(size(sol.ts));
    for ( tind = 1:length(sol.ts) )
            Lyas(tind) = (1/2)*integrate((sol.ug(:, tind) - bagSol.ug(:, tind)).^2, dz);
    end
    plot(sol.ts, Lyas, 'k-');
    grid;
    % grid minor on;
    % title('Lyapunov function');
    xlabel('t'); ylabel('\Lambda');
    ylim([0, max(Lyas(sol.ts > max(sol.ts)/5))*sqrt(2)]);
    % set(gca, 'ytick', 1.2:0.025:1.35);

end

function res = integrate(vs, dz)
    res = ( 0.5 * vs(1) + sum(vs(2:end-1)) + 0.5 * vs(end) ) * dz;
end

