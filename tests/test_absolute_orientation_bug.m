clear, close all

MWFQ=35.5; % define microwave frequency


% Spin System that defines orientation selection - SysE
clear Sys
frame1 = [0 0 0]/180*pi;
Sys.S=1/2;
Sys.g=[2.04598 2.18804];
Sys.gFrame=[frame1];
Sys.Nucs='Cu';
Sys.A=[50 612];
Sys.AFrame=[frame1];
Sys.Q = [-32];
Sys.QFrame=[frame1];
Sys.lwEndor=15;
Sys.lw=1;



%%
%Options
Opt.nKnots=10001; 
Opt.Threshold.Probe=1e-4;
Opt.Threshold.Pump=1e-4;
Opt.Output='summed';

% Define EDNMR parameters
Exp.mwFreq=MWFQ;
Exp.ExciteWidth=15;
Exp.nPoints=2048*2;
Exp.Range = [-500 0];
Exp.nu1=6*0.44; % nu1 in MHz, the scaling here comes from the fact that I used a gaussian pulse.
Exp.Tm=1;      % decay time of ENDMR nutations in us
Exp.tHTA=10;    % HTA pulse length in us
Exp.Q=80;        % Q0 of the cavity (set 1 for no frequency dependence)


Exp.Field=1176;


%%
h=figure(1);
clf
hold on

% Calculate EDNMR spectrum of
[f_sim,nmr] =  horseradish(Sys,Exp,Opt);
nmr=nmr/max(nmr);
plot(f_sim,nmr,'k','linewidth',2)


%%
frame1 = [0 50 0]/180*pi;
Sys.S=1/2;
Sys.gFrame=[frame1];
Sys.Nucs='Cu';
Sys.AFrame=[frame1];
Sys.QFrame=[frame1];
[f_sim,nmr] =  horseradish(Sys,Exp,Opt);
nmr=nmr/max(nmr);
plot(f_sim,nmr,'r','linewidth',1)
