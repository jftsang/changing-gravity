Frs = [0.1, 0.5, 1, 2];
styles = {'b', 'm', 'r', 'g'};
labels = {'Fr = 0.1', 'Fr = 0.5', 'Fr = 1', 'Fr = 2'};
nonlinearSolverConfig;

clear solFrs;
for ind = 1:length(Frs)
    u0 = @(z) Frs(ind) * (1 - 2/pi*0.05 + 0.05*sin(pi*z/2)); % almost but not quite a plug flow
    quiet = false;
    nonlinearSolverWorkhorse;
    solFrs(ind) = sol;
end

plotmSolsChis(solFrs, styles);

% figure;
% for ind = 1:length(Frs)
%     % ax = subplot(3,2,ind);
%     plotSolChi(solFrs(ind), styles{ind}, ax);
%     hold on;
% end
% hold off;
% xlim([0, 50]);
