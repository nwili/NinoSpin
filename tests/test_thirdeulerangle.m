% 14N-Nitroxide + Deuterium EDNMR spectra simulation


clear, close all

% Spin System that defines orientation selection - SysE
Sys.S=[1/2];
Sys.g =[2.04 2.08 2.19];
% Sys.D=240;
% Sys.DFrame=[0 pi/2 0];
Sys.Nucs='63Cu,63Cu';
A1=[100 300 610]/2;
Sys.A=[A1;
       A1];
Sys.Q=-35*[1 1];
% Sys.D=[240];
% Sys.DFrame=[0 0 0]/180*pi;
Sys.lw=0.03;
Sys.lwEndor=1;

Exp.mwFreq=0.3; % define microwave frequency
Exp.nPoints=2048*4;
Exp.Harmonic=0;

Exp.Range=[0 150];


angles =linspace(0,180,10);

figure(1)
clf
hold on
for i_angle = 1:numel(angles)
    
Exp.CrystalOrientation = [10 45 angles(i_angle)]/180*pi;
[x,y]=pepper(Sys,Exp);
plot(x,y)

end






%%
Exp.CrystalOrientation = [10 45 0]/180*pi;
Exp.Field=14.89;
Exp.ExciteWidth=14;

Opt.Threshold.Probe=0;
Opt.Threshold.Pump=0;


% Define EDNMR parameters
Exp.Range = [-1500 1000];
Exp.nu1=5; % nu1 in MHz
Exp.Tm=2;      % decay time of ENDMR nutations in us
Exp.tHTA=2;    % HTA pulse length in us
Exp.Q=60;        % Q0 of the cavity (set 1 for no frequency dependence)





  
figure(2)
clf
hold on
% Calculate EDNMR spectrum
sum=zeros(1,Exp.nPoints);
for i_angle = 1:numel(angles)
    
Exp.CrystalOrientation = [10 45 angles(i_angle)]/180*pi;
[x,y]=horseradish(Sys,Exp);
plot(x,y)

sum=sum+y;
end

%
plot(x,sum/numel(angles),'r','linewidth',2)


%%
angles = [0 90];
sum=zeros(1,Exp.nPoints);
for i_angle = 1:numel(angles)
    
Exp.CrystalOrientation = [10 45 angles(i_angle)]/180*pi;
[x,y]=horseradish(Sys,Exp);
% plot(x,y)

sum=sum+y;
end

plot(x,sum/numel(angles),'--b','linewidth',2)

