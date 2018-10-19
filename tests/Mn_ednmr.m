clear, close all
Sys.S=[5/2 5/2];
Sys.g=[2.00122 2.00122];
Sys.Nucs='55Mn,55Mn';
Sys.A=[-252 0 ;0 -254];
Sys.J=0;
Sys.lw=1.5;
Sys.lwEndor=1;

Sys.D = [280 220];  % D, in MHz
Sys.E=[0 0];
% Sys.DStrain=[150 150];


% Exp.mwFreq = 35;
% Exp.Range=[1100 1400];
Exp.mwFreq = 94;
Exp.Range=[3250 3500];
Exp.Harmonic=0;

% Opt.Output='separate';
Opt.Output='summed';
Opt.nKnots=[11 6];
[B_sim,Bsweep_sim]=pepper(Sys,Exp,Opt);
Bsweep_sim=Bsweep_sim/max(abs(Bsweep_sim));

plot(B_sim, Bsweep_sim)

%% setting up EDNMR simulation

Exp.ExciteWidth=15;
Exp.nPoints=2048*2;
Exp.Range = [-700 700];
Exp.nu1=10; % nu1 in MHz
Exp.Tm=4;      % decay time of ENDMR nutations in us
Exp.tHTA=30;    % HTA pulse length in us
% Exp.Q=80;        % Q0 of the cavity (set 1 for no frequency dependence)
Exp.Q=700;        


Opt.nKnots=200; %for lower symmetry
Opt.Threshold.Probe=1e-3;
Opt.Threshold.Pump=1e-2;

%% 1D

figure(1612)
clf
hold on
% Exp.Range = [0 200];
% Exp.Field=1254; 
Exp.Range = [0 200];
Exp.Field=3360; 

[freq_1d,spec_1d] =  horseradish(Sys,Exp,Opt);
plot(freq_1d,spec_1d)

% Sys.S=3/2;
% [freq_1d,spec_1d] =  horseradish(Sys,Exp,Opt);
% plot(freq_1d,spec_1d,'--')
% 
% Sys.S=1/2;
% [freq_1d,spec_1d] =  horseradish(Sys,Exp,Opt);
% plot(freq_1d,spec_1d,'--')



nu_Mn = nucgval('55Mn')*nmagn*Exp.Field*1e-3/planck*1e-6;
line1 = abs(mean(Sys.A(1,:)))/2 + nu_Mn;
line2 = abs(mean(Sys.A(1,:)))/2 - nu_Mn;
plot(line1*[1 1],ylim,line2*[1 1],ylim)
line1 = abs(mean(Sys.A(1,:))) + 2*nu_Mn;
line2 = abs(mean(Sys.A(1,:))) - 2*nu_Mn;
plot(line1*[1 1],ylim,line2*[1 1],ylim)

% nu_13C = nucgval('13C')*nmagn*Exp.Field*1e-3/planck*1e-6;
% A_C = 1.25;
% plot(nu_13C+A_C/2*[1 1],ylim); plot(nu_13C-A_C/2*[1 1],ylim);

% %% 2D
% field_pos_vec=1200:0.25:1229;
% spec_sim=zeros(numel(field_pos_vec),Exp.nPoints);
% 
% 
% tic
% for i_field=1:numel(field_pos_vec)
%    
% Exp.Field=field_pos_vec(i_field);    
%     
% % Calculate EDNMR spectrum
% 
% [freq_sim,y] =  horseradish(Sys,Exp,Opt);
% 
% i_field/numel(field_pos_vec)*100
% 
% spec_sim(i_field,:)=y;
% 
% 
% end
% sim_time = toc;
% 
% save('Mn_test')
% 
% %%
% 
% load Mn_test
% %gaussian "blur" filter, just makes the contour plot look nicer
% filt_fwhm1=10;
% filt_fwhm2=0.5;
% filt_width=50;
% x_filt=-filt_width:filt_width;
% y_filt=-filt_width:filt_width;
% [x_filt,y_filt]=meshgrid(x_filt,y_filt);
% H=normpdf(x_filt,0,filt_fwhm1).*normpdf(y_filt,0,filt_fwhm2);
% H=H/sum(H(:));
% spec_sim=filter2(H,spec_sim);
% spec_sim=spec_sim/(max(spec_sim(:)));
% 
% [freq_sim,B_ax] = meshgrid(freq_sim,field_pos_vec);
% 
% cuts = [min(spec_sim(:))*[1 0.03] [0.03 1]*-min(spec_sim(:))];
% contlevels=[linspace(cuts(1),cuts(2),20)  linspace(cuts(3),cuts(4),20) ];
% 
% 
% 
% h=figure(2);
% clf
% 
% h_bottom=subplot(4,1,[2:4]);
% hold on
% 
% [~,h_contourf]=contourf(B_ax,freq_sim,spec_sim*100,'levellist',contlevels*100','edgecolor','none');
% colormap(bluewhitered)
% %hack the colormap
% colmap=h.Colormap;
% i_white = find(sum(colmap')==3);
% colmap(-[1:2]+min(i_white),:)=[1 1 1;1 1 1];
% colormap(colmap);
% 
% 
% 
% box on
% axis([1196 1232 -750 750])
% xlimits=xlim;
% xlabel('B / mT')
% % set(gca,'XTick',[1150:30:1240])
% 
% % set(gca,'YTick',[],'YTicklabels',[])
% 
% h_col=colorbar;
% ylabel(h_col,{'Hole Depth / a.u.'})
% 
% 
% h_top=subplot(4,1, 1);
% hold on
% plot(B_sim,Bsweep_sim,'k')
% 
% 
% set(gca,'YTick',[],'XTicklabels',[],'visible','off')
% h_ax=gca; h_ax.YAxis.Visible='off'; h_ax.XAxisLocation='top';
% axis([xlimits -1.04 1.04])
% 
% 
% h=niwi_niceplot(h,[20 14]);
% 
% 
% drawnow
% h_top.Position(2)=h_top.Position(2)-0.05;
% h_top.Position(4)=h_top.Position(4)+0.08;
% drawnow
% h_top.Position(3)=h_bottom.Position(3);
% drawnow
% 
% h_bottom.Children.LineWidth=0.5;
% h_col.LineWidth=1;
% 
% set(findall(h,'-property','Fontsize'),'FontName', 'Century Gothic'); % change fontsize
% set(findall(h,'-property','Fontsize'),'Fontsize',12); % change fontsize
% 
% h.Renderer='painters';