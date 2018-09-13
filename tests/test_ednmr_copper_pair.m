% 14N-Nitroxide + Deuterium EDNMR spectra simulation


clear, close all

% Spin System that defines orientation selection - SysE
Sys.S=[1];
Sys.g =[2.04 2.19];
Sys.D=240;
Sys.DFrame=[0 pi/2 0];
Sys.Nucs='63Cu,63Cu';
A1=[100 610]/2;
Sys.A=[A1;
       A1];
Sys.Q=-35*[1 1];
% Sys.D=[240];
% Sys.DFrame=[0 0 0]/180*pi;
Sys.lw=3.9;
Sys.lwEndor=45;

Exp.mwFreq=35.5; % define microwave frequency
Exp.nPoints=2048;
Exp.Harmonic=0;

Exp.Range=[1100 1260];

pepper(Sys,Exp)
%%

% Exp.Field=1214;
Exp.ExciteWidth=14;


Opt.nKnots=51; %for lower symmetry
% Opt.nKnots=1001; %higher symmetry
Opt.Threshold.Probe=1e-3;
Opt.Threshold.Pump=1e-2;


% Define EDNMR parameters
Exp.Range = [-700 700];
Exp.nu1=5; % nu1 in MHz
Exp.Tm=2;      % decay time of ENDMR nutations in us
Exp.tHTA=2;    % HTA pulse length in us
Exp.Q=60;        % Q0 of the cavity (set 1 for no frequency dependence)


%% 2D calculation

field_pos_vec=linspace(1120,1260,28);
nmr_2D=zeros(numel(field_pos_vec),Exp.nPoints);



xsaver=cell(size(field_pos_vec));
tic
for i_field=1:numel(field_pos_vec)
   
Exp.Field=field_pos_vec(i_field);    
    
% Calculate EDNMR spectrum

[xsaver{i_field},y] =  horseradish(Sys,Exp,Opt);

i_field/numel(field_pos_vec)*100

nmr_2D(i_field,:)=y;


end
toc
% save cu_pair_test

%%

load cu_pair_test

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
xax=field_pos_vec;
[X,Y]=meshgrid(xax,yax);
colormap((jet))



% cuts = [-10 -1 1.8 7];
cuts = [min(nmr_2D(:))*[1 0.1] [0.005 0.8]*-min(nmr_2D(:))];
levels=[linspace(cuts(1),cuts(2),15)  linspace(cuts(3),cuts(4),15) ];


[~,h_contourf]=contourf(xax,yax,nmr_2D','levellist',levels);
colormap(jet(41))
colmap=h.Colormap;

whiteind=22;
mask=normpdf(1:size(colmap,1),whiteind,2)'; mask=mask/max(mask); mask=(mask*10+1); mask(whiteind,:)=100;
mask=repmat(mask,[1 3]);

colmap=1-(1-colmap)./mask;
[~,i_bg]=max(sum(colmap'));
colmap(i_bg-1:i_bg+1,:)=1;


h.Colormap=colmap;

% h_col=colorbar;
% h_col.Title.String = 'Hole Depth';
% axis([260 360 -700 700])
box on
set(gca,'Linewidth',2);

xlabel('B / mT')
ylabel('\nu_{HTA}-\nu_{obs}')

% h.Renderer='painters';

