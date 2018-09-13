clear, close all

% Spin System that defines orientation selection - SysE
Sys.S=1/2;
Sys.g=[2.04598 2.1];
Sys.Nucs='63Cu';
Sys.A=[400 612]/10;
Sys.AFrame=[0 0 0]/180*pi;
Sys.Q = -32;
Sys.lwEndor=2;
Sys.lw=5;

% Define EDNMR parameters
Exp.mwFreq=35.5;
Exp.ExciteWidth=50;
Exp.nPoints=2048;
Exp.nu1=2; % nu1 in MHz
Exp.Tm=1;      % decay time of ENDMR nutations in us
Exp.tHTA=10;    % HTA pulse length in us
Exp.Q=55;        % Q0 of the cavity (set 1 for no frequency dependence)
Exp.Temperature=300;
Exp.Harmonic=0;
% pepper(Sys,Exp)

%%
Exp.Field=1225;
Exp.Range = [-500 500];

%Options
Opt.Symmetry = symm(Sys); 
Opt.nKnots=301; 
Opt.Threshold.Probe=1e-10;
Opt.Threshold.Pump=1e-10;

%% tests
[nu,spec] = horseradish(Sys,Exp,Opt);

figure(1)
clf
plot(nu,spec)