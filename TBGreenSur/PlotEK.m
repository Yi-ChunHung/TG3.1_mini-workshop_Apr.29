%% Plotting %%
% Ef = 5.3384;
load Ek.mat

w    = 10;      % weight for the intensity of the spectrum
Name = 'title of the plot';
cax  = [-400 400]; % the range of the color map


ncolor = 1e4;
% map      = [[linspace(1,1,ncolor)' linspace(0,1,ncolor)' linspace(0,1,ncolor)'];...
map    = [linspace(1,0,ncolor)' linspace(1,0,ncolor)' linspace(1,1,ncolor)']; 
figure('position',[150 0 850 660],'paperposition',[0.25 0.25 8 10.5],...
       'papertype','usletter','numbertitle','off',...
       'PaperPositionMode','manual','paperorientation','portrait',...
       'color','w');
%pcolor(KK(2:end-1,2:end),EE(2:end-1,2:end)-Ef,...
%       ((As(2:end-1,2:end)*w))/pi),shading interp;
pcolor(KK(2:end-1,2:end),EE(2:end-1,2:end)-Ef,...
        (((As(2:end-1,2:end)-Ab(2:end-1,2:end))*w))/pi),shading interp;
colormap(map);
caxis(cax);
colorbar('Ticks',cax,'TickLabels',{'bulk','surface'});
box on

title(Name)
ylabel('\bf{Energy (eV)}','FontSize',24,'interpreter','LaTex');

ax = gca;
ax.FontSize   = 20;
ax.FontWeight = 'bold';
ax.TickLength = [0.02 0.02];
ax.XTick      = [2 Nk 2*Nk];
ax.XTickLabel = {'$\Gamma$' 'X' 'L'};
ax.LineWidth  = 2;
ax.TickLabelInterpreter='latex';
%ax.ytick([-1:0.01:1])
