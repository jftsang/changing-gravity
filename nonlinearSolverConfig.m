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

%% From DPM simulations
% mu    = @(I) tan(5.1*deg) * tanh(I/0.01)          + 0.47 * abs(I)^(0.48) * sign(I);
% dmudI = @(I) tan(5.1*deg) / 0.01 * sech(I/0.01)^2 + 0.48 * 0.47 * abs(I)^(0.48-1);

%% tanh fit
% mu    = @(I) tan(10*deg)*tanh(I/0.08) ...
%                + ( tan(20*deg) - tan(10*deg) )*tanh(I/1);
% dmudI = @(I) tan(10*deg)*sech(I/0.08)^2 / 0.08 ...
%                + ( tan(20*deg) - tan(10*deg) )*sech(I/1)^2 / 1;

%% power law
% alpha = 1/3;
% mu = @(I) tan(15*deg) * I^alpha;
% dmudI = @(I) tan(15*deg) * I^(alpha-1) * alpha;

mu = @(I) tan(15*deg) + log(I);
dmudI = @(I) 1/I;

%%%% Specify how g and theta depend on time

%% For rotation
% g = @(t) 1;
% theta = @(t) 20 * deg;
% theta = @(t)  ( 14 - 28*t/40 ) * deg; 
% theta = @(t) 16 * deg * cos(2*pi*t / 1000);
% theta = @(t) 10 * deg + (8 * deg) * cos(2*pi*(t - 200)/50) ;
% theta = @(t) 16 * deg * (1 - 0.125*erf((t-300)/8));
% g = @(t) 1;
% theta = @(t) 12.5*deg - 2.5*deg*erf((t - 30)/5);

%% For time-dependent Blasius - constant g, theta; plug flow IC
g = @(t) 1;
theta = @(t) 15*deg;


%%%% Specify initial conditions
%% Nominally a plug flow, but need a tiny bit of shear to handle singularity
u0 = @(z) 12 + 0.05*sin(0.5*pi*z);

%% Nondimensionalised grain size (equivalent to rescaling time)

% delta = 1/50;
delta = 1;

%%%% Grid size and timestepping
nz = 64;
tmax = 60;
dt = 0.1;

%%%% Options for Newton--Raphson iteration
%% Use Newton--Raphson iteration to solve implicit evolution 
%% (using Crank--Nicolson)
maxNRIter = 10;
NRPrecisionGoal = 1e-8;
NRRelaxation = 1.0;

nonlinearSolverInitialisation;
