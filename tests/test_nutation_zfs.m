%script to make sure that the Probability normalization is correct.
clear all, close all


Sys.g=gfree;
Sys.S=3/2;
Sys.D=100;
Sys.lw=1;
Sys.lwEndor=10;

Exp.mwFreq=35;
Exp.nPoints=2048;
Exp.Harmonic=0;
Exp.CrystalOrientation=[0 0 0]/pi*180;


%% EDFS
Exp.Range=[1230 1270];
% pepper(Sys,Exp)

%% EDNMR: observe CT, pump satellite

Exp.Field=1249;
Exp.ExciteWidth=10;
Exp.nu1=0.1;
Exp.Tm=10000;
Exp.Range=[-500 500];
Exp.Q=1;


nut_vec = 0:0.3:10;
spec_nut=zeros(numel(nut_vec),Exp.nPoints);

for i_nut=1:numel(nut_vec)
    Exp.tHTA=nut_vec(i_nut);
    [freq,spec_nut(i_nut,:)]= horseradish(Sys,Exp);
end

[~,poi]=min(abs(freq-200));

y=spec_nut(:,poi);


figure(1)
clf
hold on
plot(nut_vec,y);
plot(0.5/(sqrt(3)*Exp.nu1)*[1 1],ylim)
