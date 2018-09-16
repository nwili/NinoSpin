% 14N-Nitroxide + Deuterium EDNMR spectra simulation


clear, close all

% Spin System that defines orientation selection - SysE
Sys.S=1/2;
Sys.g=[2.00905506 2.00617263 2.002270418];
Sys.Nucs='14N';
Sys.A=[4.2 4.2 34.1]*2.8;
Sys.Q=[1.2 0.5 -1.7];
Sys.lw = 3;
Sys.lwEndor=1.4;

Exp.Field=1220;
Exp.nPoints=2048;
Exp.Harmonic=0;


pepper(Sys,Exp)
%%

% Exp.Field=1214;
Exp.ExciteWidth=4;


Opt.nKnots=51; %for lower symmetry
% Opt.nKnots=1001; %higher symmetry
Opt.Threshold.Probe=1e-4;
Opt.Threshold.Pump=1e-4;


% Define EDNMR parameters
Exp.Range = [-70 70];
Exp.nu1=2; % nu1 in MHz
Exp.Tm=4;      % decay time of ENDMR nutations in us
Exp.tHTA=2;    % HTA pulse length in us
Exp.Q=60;        % Q0 of the cavity (set 1 for no frequency dependence)

%% 1D trial
Exp.mwFreq=34.15;

[nu,spec] = cheesy(Sys,Exp,Opt);

plot(nu,spec)



%% 2D calculation

mw_pos_vec=linspace(34.09,34.33,60);
nmr_2D=zeros(numel(mw_pos_vec),Exp.nPoints);



xsaver=cell(size(mw_pos_vec));
tic
for i_field=1:numel(mw_pos_vec)
   
Exp.mwFreq=mw_pos_vec(i_field);    
    
% Calculate EDNMR spectrum

[xsaver{i_field},y] =  cheesy(Sys,Exp,Opt);

i_field/numel(mw_pos_vec)*100

nmr_2D(i_field,:)=y;


end
toc
save cheesy_test

%%

load cheesy_test

% for ii=1:size(nmr_2D)
%    nmr_2D(ii,:)=nmr_2D(ii,:)/max(nmr_2D(ii,:)); 
% end


% gaussian "blur" filter, just makes the contour plot look nicer
filt_fwhm1=10;
filt_fwhm2=1;
filt_width=50;
x_filt=-filt_width:filt_width;
y_filt=-filt_width:filt_width;
[x_filt,y_filt]=meshgrid(x_filt,y_filt);
H=normpdf(x_filt,0,filt_fwhm1).*normpdf(y_filt,0,filt_fwhm2);
H=H/sum(H(:));
spec_proc=filter2(H,nmr_2D);



nmr_2D=spec_proc/max(spec_proc(:))*100;
% nmr_2D=nmr_2D-median(nmr_2D(:));


h=figure(2);
clf
hold on
% box on
yax=xsaver{1};
xax=mw_pos_vec;
[X,Y]=meshgrid(xax,yax);
colormap((jet))



% cuts = [-10 -1 1.8 7];

levels=[linspace(2,100,15)   ];


[~,h_contour]=contour(xax,yax,nmr_2D','levellist',levels);
colormap(jet(41))
colmap=h.Colormap;


