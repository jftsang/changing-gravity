pkg load statistics symbolic;

nonlinearSolverConfig;
tic;

%% Find the Bagnoldian solution
bagSol = nonlinearSolverBagnoldians(mu, g, theta, ts, zs);
toc;

%% Go!
nonlinearSolverWorkhorse;
nonlinearSolverDerivatives;
toc;

%% Plot the solution and compare it to the Bagnoldian solution
plotSol(sol);
toc;
