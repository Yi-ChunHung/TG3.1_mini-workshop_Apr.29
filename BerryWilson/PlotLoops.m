%%% Choose the kz in linspace(-0.5,0.5,ibz)
ibz = 1;

%%
load chern.mat

b3(ibz)

figure('Name','Chern')
hold on
title('Berry Phase')
plot(b2,fff(:,ibz)','b.','MarkerSize',15)
axis([b2(1),b2(end),-pi,pi]);
set(gca,'ytick',[-pi,0,pi],'yticklabel',{'-\pi','0','\pi'},'fontsize',20);
if mode
    set(gca,'xtick',[b2(1),0,b2(end)],'xticklabel',{'-\pi','0','\pi'},'fontsize',20);
else
    set(gca,'xtick',[0,b2(end)],'xticklabel',{'0','\pi'},'fontsize',20);
end
set(gca,'linewidth',2,'FontSize',20)
ylabel('\Theta')
xlabel('k_y')
box on
axis square

%% distinguish the Chern number
Chern_temp   = unwrap(fff(:,ibz));
Chern_number = (Chern_temp(end) - Chern_temp(1))/(2*pi);
fprintf('Chern number = %f \n',Chern_number)


%%
load Z2.mat

figure('Name','Z2')
hold on
title('Hybrid Wannier Centers')
% plot(by,theta_Dm,'b.','MarkerSize',15)
plot(b2,theta_Dm(:,:,ibz),'k.','MarkerSize',15);
axis([b2(1),b2(end),-pi,pi]);
set(gca,'ytick',[-pi,0,pi],'yticklabel',{'-\pi','0','\pi'},'fontsize',20);
if mode
    set(gca,'xtick',[b2(1),0,b2(end)],'xticklabel',{'-\pi','0','\pi'},'fontsize',20);
else
    set(gca,'xtick',[0,b2(end)],'xticklabel',{'0','\pi'},'fontsize',20);
end
set(gca,'linewidth',2,'FontSize',20)
ylabel('\theta')
xlabel('k_y')
box on
axis square
