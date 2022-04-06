load Ek.mat

Name = {'010'};
cax  = [-100 100];
AA   = As - Ab; % Ab:bulk As:surface
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
map    = [[linspace(1,1,ncolor)' linspace(0,1,ncolor)' linspace(0,1,ncolor)'];...
          [linspace(1,0,ncolor)' linspace(1,0,ncolor)' linspace(1,1,ncolor)']]; 
colormap(map);
caxis(cax);
colorbar('Ticks',cax,'TickLabels',{'bulk','surface'});
% colorbar('Ticks',cax,'TickLabels',{'dn','up'});

pcolor(kd,wef,AA/pi*100),shading interp

axis square,box on ,set(gca,'linewidth',2,'FontSize',16)
plot([kd(1) kd(end)],[0 0],'k-.')
ax = gca;
ax.XTick      = [0 nk/2 nk];
ax.XTickLabel = {'k_x=-\pi' 'k_x=0' 'k_x=\pi'};
%ax.XTickLabel = {'k_x=0.505\pi' 'k_x=0.5075\pi' 'k_x=0.51\pi'};

title(Name)

axis([0 kd(end) ymin ymax])
ylabel('Energy (eV)')