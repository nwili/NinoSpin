clear all, close all


%close up
file = '20180803_1059_fstep_EDNMR';

% options for uwb_eval
clear options;
options.plot = 0; % 0 for do not plot by uwb_eval
options.evlen = 1000; % enforce a 256 point evaluation window
options.cor_phase = 0; %0; %1; %1; % phase all data points, since amplitude sweeps make no sense without phasing

% get downconverted data by uwb_eval
output = uwb_eval(file,options);

expLO=output.exp.LO;

dta_cont = output.dta_avg; % averaged echo transients
dta_ev = output.dta_ev; % echo integrals
ev_coll = dta_ev;
dta_x_cont = output.dta_x{1}; % axis of first indirect dimension
t_ax = output.t_ax;
exp_cont = output.exp;

cent_frq=expLO+1.5;

x=dta_x_cont+expLO;
[x,ind]=sort(x);
y=real(ev_coll);
y=y(ind);

% y=datasmooth(y,2,'savgol',1,0);
y=datasmooth(y,3,'binom');
% y = medfilt1(y,20);
y=1-(y/max((y)));
spec_exp=y/max(y);


f_exp =1e3*(x-cent_frq);

%%
% 1Copper EDNMR spectra simulation
MWFQ=35.5; % define microwave frequency

% Spin System that defines orientation selection - SysE
Sys.S=1/2;
Sys.g=[2.04598 2.18804];
% Sys.gStrain=[0.0064 0.0038];
Sys.Nucs='63Cu';
Sys.A=[50 612];
Sys.AFrame=[0 0 0]/180*pi;
% Sys.HStrain = [100 100];
Sys.Q = -32;
Sys.lwEndor=0;

%Options
Opt.Symmetry = symm(Sys); 
Opt.nKnots=10001; 
Opt.Threshold.Probe=1e-4;
Opt.Threshold.Pump=1e-4;

% Define EDNMR parameters
Exp.mwFreq=MWFQ;
Exp.ExciteWidth=20;
Exp.nPoints=2048*2;
Exp.Range = [-500 0];
Exp.nu1=2; % nu1 in MHz
Exp.Tm=1.5;      % decay time of ENDMR nutations in us
Exp.tHTA=10;    % HTA pulse length in us
Exp.Q=55;        % Q0 of the cavity (set 1 for no frequency dependence)


% load resprofile_sca10000;
% pulse_scale = 0.02;
% Exp.ResProfile.freq = nut_y;
% Exp.ResProfile.w1 = nu1*pulse_scale*1e9;

Exp.Field=1176;

% Calculate EDNMR spectrum of Cu63
[f_sim,y] =  horseradish(Sys,Exp,Opt);
y=y/max(y);
nmr_63=y;

% Calculate EDNMR spectrum of Cu63
Sys.Nucs='65Cu';
Sys.A=Sys.A*nucgval('65Cu')/nucgval('63Cu');
Sys.Q=Sys.Q*nucqmom('65Cu')/nucqmom('63Cu');


[~,y] =  horseradish(Sys,Exp,Opt);
y=y/max(y);
nmr_65=y;


%%
figure(1)
clf
hold on

ch = 100:950;
f_exp=f_exp(ch);
spec_exp=spec_exp(ch);

% plot(f_sim,nmr_65*nucabund('65Cu'))
% plot(f_sim,nucabund('63Cu')*nmr_63)

spec_sum_sim=nmr_65*nucabund('65Cu')+nucabund('63Cu')*nmr_63;
spec_sum_sim=spec_sum_sim/max(spec_sum_sim);

plot(f_sim,spec_sum_sim/2,'r','linewidth',2)

[spec_sum_sim,lw_min,alpha_min] = rescale_n_lw(f_sim,spec_sum_sim,f_exp,spec_exp,'lsq1');
% spec_exp = rescale_n_lw(spec_exp,spec_sum_sim,'lsq1');


plot(f_exp,spec_sum_sim,'r','linewidth',2)
plot(f_exp,spec_exp,'k')
plot(f_exp,3*(spec_sum_sim-spec_exp')-1)

lw_min
alpha_min

% axis([-400 -120 ylim])