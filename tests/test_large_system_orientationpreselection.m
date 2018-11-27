clear, close all

Exp.mwFreq=35.5; % define microwave frequency
Exp.CrystalOrientation=[0 10 0;
                        0 44 0;]/180*pi;
Exp.CrystalWeight=[100 1];
Exp.Harmonic=0;

% Spin System that defines orientation selection - SysE
Sys.S=1/2;
Sys.g=[2.04598 2.18804];
% Sys.gStrain=[0.0064 0.0038];
Sys.Nucs='63Cu';
Sys.A=-[50 612];
Sys.AFrame=[0 0 0]/180*pi;
% Sys.HStrain = [100 100];
Sys.Q = 32;
Sys.lw=1;
Sys.lwEndor=10;


pepper(Sys,Exp)
%%

%Options
Opt.Symmetry = symm(Sys); 
Opt.nKnots=1001; 
Opt.Threshold.Probe=1e-4;
Opt.Threshold.Pump=1e-4;

% Define EDNMR parameters
Exp.ExciteWidth=20;
Exp.nPoints=2048*2;
Exp.Range = [-500 500];
Exp.nu1=5; % nu1 in MHz
Exp.Tm=1.5;      % decay time of ENDMR nutations in us
Exp.tHTA=10;    % HTA pulse length in us
Exp.Q=55;        % Q0 of the cavity (set 1 for no frequency dependence)




Exp.Field=1189;
% Calculate EDNMR spectrum
[f_sim,y] =  horseradish(Sys,Exp,Opt);


figure(2)
clf
hold on
plot(f_sim,y)