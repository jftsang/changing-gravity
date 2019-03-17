%%%% Initialisation

%% Grid points
zs = linspace(0, 1, nz); 
dz = zs(2);
dzs = zs(2:nz) - zs(1:nz-1); dzs = [dzs, 0];

%% Half-grid points
zphalf = zs + dz/2;
zmhalf = zs - dz/2;

%% Time 
ts = 0:dt:tmax;
nt = length(ts);

[tg, zg] = meshgrid(ts, zs);

%% Mesh grid for the solution
ug = zeros(nz, nt);

%% Initial conditions
ug(:,1) = arrayfun(u0, zs);
