pkg load statistics symbolic;

% mu = @(I) tan(10*pi/180) * tanh(I/0.02) + I;
% dmudI = @(I) tan(10*pi/180) * sech(I/0.02)^2 / 0.02 + 1;

%% Real particles
%{
mu = @(I) 0.0746 + 0.5029 * I^0.4336;
dmudI = @(I) 0.5029 * 0.4336 * I^(0.4336 - 1);
%}
%{
mu    = @(I) 0.0746 * tanh(I/0.01) + sign(I) * 0.5029 * abs(I)^0.4336;
dmudI = @(I) 0.0746 / 0.01 * sech(I/0.01)^2 + 0.4336 * 0.5029 * abs(I)^(0.4336 - 1);
%}
% Not working well --- have infinite dmu/dI at I = 0
% Might be actually well-posed, but numerically no good at all.
%{
mu    = @(I) sign(I) * 0.5029 * abs(I)^0.4336;
dmudI = @(I) 0.4336 * 0.5029 * abs(I)^(0.4336 - 1);
%}

%% Mollified affine fit
% mu    = @(I) tan(10 * pi/180) * tanh(I/0.03) + 0.3585*I;
% dmudI = @(I) tan(10 * pi/180) / 0.03 * sech(I/0.03)^2 + 0.3585;

%% From DPM simulations
mu    = @(I) tan(5.1*pi/180) * tanh(I/0.01)          + 0.47 * abs(I)^(0.48) * sign(I);
dmudI = @(I) tan(5.1*pi/180) / 0.01 * sech(I/0.01)^2 + 0.48 * 0.47 * abs(I)^(0.48-1);

% Simple power law
% mu = @(I) I;
% dmudI = @(I) 1;

%%%% g and theta

%% For rotation
% g = @(t) 1;
% theta = @(t) 20 * pi/180;
% theta = @(t)  ( 14 - 28*t/40 ) * pi/180; 
% theta = @(t) 16 * pi/180 * cos(2*pi*t / 1000);
% theta = @(t) 16 * pi/180 * (1 - 0.125*erf((t-300)/8));
% theta = @(t) 10 * pi/180 + (8 * pi/180) * cos(2*pi*(t - 200)/50) ;

%% For time-dependent Blasius - constant g, theta; plug flow IC
g = @(t) 1;
theta = @(t) 16 * pi/180;
%% Quasi-plug flow, but need a tiny bit of shear to handle singularity
u0 = @(z) 1.00 + 0.005*sin(0.5*pi*z);
% u0 = @(z) z;

% u0 = @(z) 0.5 + 0.5 * tanh((z-0.5)/0.1);
% u0 = @(z) 0.19 * ( 1 - (1-z)^(3/2) );
% u0 = @(z) 0.2137 * ( 1 - (1-z)^(3/2) );
% u0 = @(z) 0.15 * (tanh((z-0.333)/0.1) - tanh((z-0.667)/0.1));


% delta = 1/50;
delta = 1;

nz = 48;
tmax = 20;
dt = 0.01;

% Use Newton--Raphson iteration to solve implicit evolution 
% (using Crank--Nicolson)
maxNRIter = 10;
NRPrecisionGoal = 1e-8;
NRRelaxation = 1.0;

nonlinearSolverInitialisation;
