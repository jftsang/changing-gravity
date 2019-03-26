nonlinearSolverConfig;

%% Initialisation
zs = linspace(0, 1, nz); 
dz = zs(2);
% Half-grid points
zphalf = zs + dz/2;
zmhalf = zs - dz/2;

nNRIters = 0;
NRChange = inf;

usteady = arrayfun(u0, zs)';

t = 0;
gnow = g(t);
thetanow = theta(t);

tic;
while (true)
    %% Vector of residuals and Jacobian matrix
    Res = zeros(nz, 1);
    Jac = zeros(size(nz));

    %% No-slip condition
    Res(1) = usteady(1);
    Jac(1,1) = 1;

    for j = 2:nz-1
      % Inertial number at half-points
      Ip = (usteady(j+1) - usteady(j)) ...
             / (dz * ( gnow * cos(thetanow) * (1-zphalf(j)) )^(1/2) );
      Im = (usteady(j) - usteady(j-1)) ...
             / (dz * ( gnow * cos(thetanow) * (1-zmhalf(j)) )^(1/2) );

      % Derivative of inertial number at half-points

             
      Res(j) =  ...
                - gnow * sin(thetanow) ...
                - gnow * cos(thetanow) * ...
                      (   mu(Ip) * (1-zphalf(j)) ...
                        - mu(Im) * (1-zmhalf(j)) ...
                      ) / dz;
              
      Jac(j,j) = (gnow * cos(thetanow))^(1/2) / (dz^2) ...
                    * ( dmudI(Ip) * (1-zphalf(j))^(1/2) ...
                      + dmudI(Im) * (1-zmhalf(j))^(1/2) );
      Jac(j,j+1) = - ( gnow*cos(thetanow))^(1/2) / (dz^2) ...
                         * dmudI(Ip) * (1-zphalf(j))^(1/2);
      Jac(j,j-1) = - ( gnow*cos(thetanow))^(1/2) / (dz^2) ...
                         * dmudI(Im) * (1-zmhalf(j))^(1/2);                     

    end
    clear Ip Im;

    %% No-stress condition
    Res(nz) = ((1/2)*usteady(nz-2) - 2*usteady(nz-1) + (3/2)*usteady(nz))/dz ;
    Jac(nz,nz-2) = 1/(2*dz);
    Jac(nz,nz-1) = -2/dz;
    Jac(nz,nz)   = 3/(2*dz);

    Jac = sparse(Jac);

    % NR iteration
    dusteady = -Jac \ Res;
    nNRIters++;
    usteady += dusteady;

    % Check if the change dusteady falls below the NRPrecisionGoal
    NRChange = max(abs(dusteady));
    if (NRChange < NRPrecisionGoal)
      break;
    end

    % Check if we have exceeded the maximum number of NR iterations  
    if (nNRIters >= maxNRIter)
      warning(sprintf( ...
        'NR failed to achieve NRPrecisionGoal %e after %d steps.', ...
        NRPrecisionGoal, nNRIters));
      break;
    end
    
end

usteadyz = zeros(nz, 1);
usteadyz(1) =  ( -(3/2)*usteady(1)  + 2*usteady(2) - (1/2)*usteady(3) )/ dz;
for j = 2:nz-1
    usteadyz(j) = ( usteady(j+1) - usteady(j-1) ) / (2*dz);
end
usteadyz(nz) = ( (3/2)*usteady(nz) - 2*usteady(nz-1) + (1/2)*usteady(nz-2) ) / dz;
isteady = usteadyz ./ (gnow * cos(thetanow) * (1 - zs') ).^(1/2);

toc;

%% Bagnoldian solution - need to invert mu(I)
Ith = fzero(@(I) mu(I) - tan(thetanow), 0);

usteadyBag = 2/3 * Ith * (gnow * cos(thetanow))^(1/2) ...
                      * ( 1 - (1-zs).^(3/2) );

%% Plots

figure;
subplot(1,3,1);
plot( usteady, zs, ...
      usteadyBag, zs );
subplot(1,3,2);
plot( usteadyz, zs, ...
      Ith * (gnow * cos(thetanow))^(1/2) * (1 - zs).^(1/2), zs);
subplot(1,3,3);
plot( isteady, zs, ...
      Ith  + zs*0, zs );
