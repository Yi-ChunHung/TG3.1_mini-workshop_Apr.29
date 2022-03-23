function PlotDispersion(Ef,Name)

load Ekplot.mat

fprintf('Start plotting... \n')
tic
%=======================================================================================================
%% --- Plot the Band Structure and the Partial Charage Distributions ---- %%%
%% --- Also, the color map of the Spin Polarization along a direction --- %%%
% ------------------------------------------------------------------------- %
hold on
for j = 1:norb
	plot(p,Ek(:,j)-Ef,'k.','LineWidth',1);
	if length(PCD_orb_1) ~= 0
		scatter(p,Ek(:,j)-Ef,100*PCD_1(:,j)+1.0e-14,'r','filled');
		%scatter(p,Ek(:,j)-Ef,100*PCD_1(:,j)+1.0e-14,PCD_1(:,j),'filled');
	end
	if length(PCD_orb_2) ~= 0
		scatter(p,Ek(:,j)-Ef,100*PCD_2(:,j)+1.0e-14,'b','filled');
		%scatter(p,Ek(:,j)-Ef,100*PCD_2(:,j)+1.0e-14,PCD_2(:,j),'filled');
	end
    if length(PCD_orb_3) ~= 0
            scatter(p,Ek(:,j)-Ef,100*PCD_3(:,j)+1.0e-14,'m','filled');
			%scatter(p,Ek(:,j)-Ef,100*PCD_3(:,j)+1.0e-14,PCD_3(:,j),'filled');
    end
    if length(PCD_orb_4) ~= 0
            scatter(p,Ek(:,j)-Ef,100*PCD_4(:,j)+1.0e-14,'g','filled');
			%scatter(p,Ek(:,j)-Ef,100*PCD_4(:,j)+1.0e-14,PCD_4(:,j),'filled');
    end
end
%colorbar % to see the quantitative partial charge distribution
hold off
% ------------------------------------------------------------------------- %
title(Name)
ylabel('\bf{Energy (eV)}','FontSize',24,'interpreter','LaTex');
%xlabel('$\vec{k}/2\pi$ in the reciprocal lattice representation','Interpreter','latex','FontSize',15);
axis([0 p(end) min(Ek-Ef,[],'all') max(Ek-Ef,[],'all')]);
%axis([0 p(end) -1 2]);
ax = gca;
ax.FontSize   = 24;
ax.FontWeight = 'bold';
ax.TickLength = [0.02 0.02];
ax.XTickLabel = symlb;
ax.XTick      = sympt;
ax.LineWidth  = 0.5;
ax.TickLabelInterpreter='latex';
grid on
set(gcf, 'Position',  [150, 150, 2000, 1600])
toc
%=======================================================================================================