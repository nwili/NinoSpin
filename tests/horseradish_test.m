clear, close all

% Spin System that defines orientation selection - SysE
Sys.S=1/2;
Sys.g=[2.04598 2.18];
Sys.Nucs='Cu,2H';
Sys.A=[80 600; [1 1.5]*0.2 ];
% Sys.AFrame=[0 0 0]/180*pi;
Sys.Q = [-32; 0];
Sys.lwEndor=20;
Sys.lw=5;

% Define EDNMR parameters
Exp.mwFreq=35.5;
Exp.ExciteWidth=14;
Exp.nPoints=2048;
Exp.nu1=10; % nu1 in MHz
Exp.Tm=1;      % decay time of ENDMR nutations in us
Exp.tHTA=10;    % HTA pulse length in us
Exp.Q=55;        % Q0 of the cavity (set 1 for no frequency dependence)
Exp.Temperature=300;
Exp.Harmonic=0;
% Exp.CrystalOrientation=[0 45 0]/180*pi;

Opt=struct;

figure(1)
clf
subplot(2,1,1)
Exp.Range=[1110 1280];
[x,y]=pepper(Sys,Exp,Opt);
plot(x,y)
hold on

%%
Exp.Field=1240;
plot(Exp.Field*[1 1],ylim)


Exp.Range = [-500 500];

%Options
Opt.nKnots=51; 
Opt.Output='separate';
% Opt.Threshold.Probe=1e-10;
% Opt.Threshold.Pump=1e-10;

%% tests
[nu,spec] = horseradish(Sys,Exp,Opt);

figure(1)
subplot(2,1,2)
plot(nu,spec)