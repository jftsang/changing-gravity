%% function fig = plotSolChiPrintable(sol)

function fig = plotSolChiPrintable(sol)
  fig = figure( ...
            'paperunits', 'centimeters', ...
            'papersize', [18 6] ...
        );
            % 'outerposition', [0 0 0 0] ...
            % 'outerposition', [0 0 0 0] ...
    % fig = figure;
    dz = sol.zs(2) - sol.zs(1);
    chis = zeros(size(sol.ts));
    for ( tind = 1:length(sol.ts) )
        chis(tind) = integrate(sol.ug(:, tind).^2, dz)  ...
                        / integrate(sol.ug(:, tind), dz)^2; 
    end
    plot(sol.ts.^(1/(1+0.48)), chis, 'k-');
    grid;
    % title('Shape factor');
    xlabel('t^{1/(1+\alpha)}'); ylabel('\chi');
    ylim([1.0, 1.35]);
    set(gca, 'ytick', 1.0:0.025:1.35);
end
