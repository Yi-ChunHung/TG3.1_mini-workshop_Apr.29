load Ek_norashba_HHD.mat
As_ml = As;
Ab_ml = Ab;

load Ek_rashba_HHD.mat

Name = {'010 overlape'};
%Name = {'010 no rashba'};
cax  = [-100 100];
% AA   = As - Ab; % Ab:bulk As:surface
% AA   = As_sz;

wef  = w - ef;
kd   = linspace(0,nk,nk);
nk   = length(kd);
%xlabel('K_x')
ymin = min(w);
ymax = max(w);
figure; 
hold on
set(gca,'XTick',[],'XTickLabel',{''})
ncolor = 1000;
map1    = [[linspace(1,1,ncolor)' linspace(0,1,ncolor)' linspace(0,1,ncolor)'];...
           [linspace(1,0,ncolor)' linspace(1,0,ncolor)' linspace(1,1,ncolor)']]; 
map2    = [[linspace(0.8,1,ncolor)' linspace(0.8,1,ncolor)' linspace(0.8,1,ncolor)'];...
           [linspace(1,0,ncolor)' linspace(1,0,ncolor)' linspace(1,0,ncolor)']];
map3    = [[linspace(0,1,ncolor)' linspace(0,1,ncolor)' linspace(0,1,ncolor)'];...
           [linspace(1,0,ncolor)' linspace(1,0,ncolor)' linspace(1,1,ncolor)']]; 
caxis(cax);
%colorbar('Ticks',cax,'TickLabels',{'bulk','surface'});
%colorbar('Ticks',cax,'TickLabels',{'dn','up'});
%colorbar('Ticks',cax,'TickLabels',{'no rashba','rashba'});

%h1 = axes();
pcolor(kd,wef,(2*As_ml - Ab_ml)/pi*100);,shading interp
colormap(map2);
freezeColors

h2 = pcolor(kd,wef,(As + 0*Ab)/pi*100);,shading interp
h2.FaceAlpha = 0.65;
%colormap(map1);
%colormap(map2);
colormap(map3);

axis square,box on ,set(gca,'linewidth',2,'FontSize',16)
plot([kd(1) kd(end)],[0 0],'k:')
ax = gca;
ax.XTick      = [0 nk/2 nk];
%ax.XTickLabel = {'k_x=-\pi' 'k_x=0' 'k_x=\pi'};
ax.XTickLabel = {'\Gamma' 'X' '\Gamma'};

%title(Name)

%axis([0 kd(end) ymin ymax])
%axis([0 kd(end) -0.1 0.1])
axis([0 kd(end) -0.5 0.5])
ylabel('Energy (eV)')

saveas(gcf,'010_overlape.png')