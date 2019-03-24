%%%% Specify the function mu(I)
deg = pi/180;

% mu = @(I) tan(10*deg) * tanh(I/0.02) + I;
% dmudI = @(I) tan(10*deg) * sech(I/0.02)^2 / 0.02 + 1;

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
% mu    = @(I) tan(10 * deg) * tanh(I/0.03) + 0.3585*I;
% dmudI = @(I) tan(10 * deg) / 0.03 * sech(I/0.03)^2 + 0.3585;

%% From DPM simulations
% mu    = @(I) tan(5.1*deg) * tanh(I/0.01)          + 0.47 * abs(I)^(0.48) * sign(I);
% dmudI = @(I) tan(5.1*deg) / 0.01 * sech(I/0.01)^2 + 0.48 * 0.47 * abs(I)^(0.48-1);

%% Jop fit
% mu    = @(I) tan(10*deg)*tanh(I/0.08) ...
%                + ( tan(20*deg) - tan(10*deg) )*tanh(I/1);
% dmudI = @(I) tan(10*deg)*sech(I/0.08)^2 / 0.08 ...
%                + ( tan(20*deg) - tan(10*deg) )*sech(I/1)^2 / 1;

%% power law
I15 = ( (tan(20*deg)-tan(10*deg)) / (tan(15*deg)-tan(10*deg)) )^(1/7);
mu = @(I) tan(10*deg) + (tan(15*deg)-tan(10*deg)) * (I/I15)^(1/7);
dmudI = @(I) (tan(15*deg)-tan(10*deg)) * (I/I15)^(-6/7) * (1/7) / I15;

%%%% Specify how g and theta depend on time

%% For rotation
% g = @(t) 1;
% theta = @(t) 20 * deg;
% theta = @(t)  ( 14 - 28*t/40 ) * deg; 
% theta = @(t) 16 * deg * cos(2*pi*t / 1000);
% theta = @(t) 16 * deg * (1 - 0.125*erf((t-300)/8));
% theta = @(t) 10 * deg + (8 * deg) * cos(2*pi*(t - 200)/50) ;

%%%% Specify initial conditions

%% For time-dependent Blasius - constant g, theta; plug flow IC
g = @(t) 1;
theta = @(t) 15*deg;
%% Quasi-plug flow, but need a tiny bit of shear to handle singularity
u0 = @(z) 2.00 + 0.05*sin(0.5*pi*z);
% u0 = @(z) z;

% u0 = @(z) 0.5 + 0.5 * tanh((z-0.5)/0.1);
% u0 = @(z) 0.19 * ( 1 - (1-z)^(3/2) );
% u0 = @(z) 0.2137 * ( 1 - (1-z)^(3/2) );
% u0 = @(z) 0.15 * (tanh((z-0.333)/0.1) - tanh((z-0.667)/0.1));

%% Nondimensionalised grain size (equivalent to rescaling time)

% delta = 1/50;
delta = 1;

%%%% Grid size and timestepping
nz = 96;
tmax = 80;
dt = 0.05;

%%%% Options for Newton--Raphson iteration
%% Use Newton--Raphson iteration to solve implicit evolution 
%% (using Crank--Nicolson)
maxNRIter = 10;
NRPrecisionGoal = 1e-8;
NRRelaxation = 1.0;

nonlinearSolverInitialisation;
