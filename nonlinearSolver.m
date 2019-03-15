nonlinearSolverConfig;
tic;

%% Find Bagnoldian solution
bagSol = nonlinearSolverBagnoldians(mu, g, theta, ts, zs);
toc;

%% Go!
nonlinearSolverWorkhorse;
nonlinearSolverDerivatives;
toc;

%% Diagram
plotSol(sol, bagSol);
toc;
