if (~exist('ts', 'var') || ~exist('ug', 'var'))
    nonlinearSolverInitialisation;
end

%% Pressure (probably lithostatic, could also be constant for a shear cell)
function p = pr(z, g, theta)
  p = g * cos(theta) * ( 1 - z );
end

%% function [Ihs, Ips, Ims] = inertials(us, zs, g, theta)
%% Calculates the inertial number at the interior gridpoints 
%% and at the half-grid points.
function [Ihs, Ips, Ims] = inertials(us, zs, g, theta)
      dz = zs(2); nz = length(zs);
      
      % Force zs to be a column vector.
      zs = zs.';


      % Half-grid points
      zphalf = zs + dz/2;
      zmhalf = zs - dz/2;
      Ihs = zeros(nz, 1); Ips = Ihs; Ims = Ihs;

      % Inertial number here
      parfor j = 2:nz-1
          % Inertial number at grid points ('here')
          Ihs(j) = (us(j+1) - us(j-1)) ...
            ./ (2 * dz * pr(zs(j), g, theta)^(1/2) );
          % Inertial number at half-points
          Ips(j) = (us(j+1) - us(j)) ...
            ./ (dz * pr(zphalf(j), g, theta)^(1/2) );
          Ims(j) = (us(j) - us(j-1)) ...
            ./ (dz * pr(zmhalf(j), g, theta)^(1/2) );
      end
end

%% Main loop
NRFailedSteps = zeros(nt, 1);      % Number of steps at which NR fails to converge
for tind = 1:(nt-1)
  % Evaluate time, g and theta at the half-time between 
  % this timestep and the next
  thalf = (ts(tind) + ts(tind+1))/2;
  ghalf = g(thalf);
  thetahalf = theta(thalf);

  nNRIters = 0;
  NRChange = inf;
  ucurr = ug(:, tind); % u at the current timestep
  utest = ug(:, tind); % a trial for the next timestep
                       %   which we improve using NR

  % Loop for NR iteration
  while (true)
%% Vector of residuals and Jacobian matrix
    Res = zeros(nz, 1);
    Jac = zeros(size(nz));

%% Basal condition - no-slip or generalised Knudsen
    % Res(1) = utest(1) - 0.20 * heaviside(500 - thalf) ;
    Res(1) = utest(1);
    Jac(1,1) = 1;
    
%% Residues at the interior points
    [Ihs, Ips, Ims] = inertials(utest, zs, ghalf, thetahalf);
    parfor j = 2:nz-1
        % Inertial number here
        Ih = Ihs(j); 
        % Inertial number at half-points (p, m == plus, minus)
        Ip = Ips(j); Im = Ims(j);

        % Derivatives of mu by finite differences
        Res(j) = (utest(j) - ucurr(j))/dt ...
                - delta * ghalf * sin(thetahalf) ...
                - delta * ( ...
                      mu(Ip) * pr(zphalf(j), ghalf, thetahalf) ...
                    - mu(Im) * pr(zmhalf(j), ghalf, thetahalf) ...
                  ) / dz;
              
        % Need dmudI for the Jacobian
        Jac(j,j) = 1/dt + delta / (dz^2) ...
                    * ( dmudI(Ip) * pr(zphalf(j), ghalf, thetahalf)^(1/2) ...
                      + dmudI(Im) * pr(zmhalf(j), ghalf, thetahalf)^(1/2) );
        Jac(j,j+1) = - delta / (dz^2) ...
                         * dmudI(Ip) * pr(zphalf(j), ghalf, thetahalf)^(1/2);
        Jac(j,j-1) = - delta / (dz^2) ...
                         * dmudI(Im) * pr(zmhalf(j), ghalf, thetahalf)^(1/2);                     
    end
    clear Ip Im;

%% No-stress condition at endpoint
    Res(nz) = ((1/2)*utest(nz-2) - 2*utest(nz-1) + (3/2)*utest(nz))/dz ;
    Jac(nz,nz-2) = 1/(2*dz);
    Jac(nz,nz-1) = -2/dz;
    Jac(nz,nz)   = 3/(2*dz);

%% No-slip condition at endpoint
%{
    Res(nz) = utest(nz) - 1;
    Jac(nz,nz) = 1;
%}
    
%% NR iteration
    Jac = sparse(Jac);
    dutest = -NRRelaxation * Jac \ Res;
    nNRIters++;
    utest += dutest;

    % Check if the change dutest falls below the NRPrecisionGoal
    NRChange = max(abs(dutest));
    %{
    printf('t = %f, dutest = %e, NRChange = %e, after %d iterations\n', ...
        thalf, dutest(nz), NRChange, nNRIters);
    %}
    if (NRChange < NRPrecisionGoal)
        ug(:, tind+1) = utest; % Fill in the next timestep. 
                             % We will use this as ucurr 
                             %   in the next loop iteration.
        clear ucurr utest;
        break;
    end

    % Check if we have exceeded the maximum number of NR iterations  
    if (nNRIters >= maxNRIter)
        NRFailedSteps(tind) = 1;  
        warning(sprintf( ...
          'At t = %f: NR failed to achieve NRPrecisionGoal %e after %d steps.\nNRChange = %e, cond(Jac) = %e', ...
          thalf, NRPrecisionGoal, nNRIters, NRChange, cond(Jac)));
        ug(:, tind+1) = utest;
        clear ucurr utest;
        break;
    end
    
  end
  printf('Completed %d of %d timesteps. NRChange = %.3e in %d iterations, cond = %e\n', ...
    tind, nt-1, NRChange, nNRIters, cond(Jac));
end

printf('NR failed to converge on %d out of %d timesteps\n', 
  sum(NRFailedSteps), nt-1);

sol = struct('ts', ts, 'zs', zs, 'tg', tg, 'zg', zg, 'ug', ug);
