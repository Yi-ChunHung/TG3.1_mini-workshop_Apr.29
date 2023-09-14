%==========================================================================
%              -- Band Structure Plot for ftn58sparse --
%
%              Hsin Lin, Tay-Rong Chang, Mohammad Yahyavi 
%==========================================================================
clear; close all; clc;   
%--------------------------------------------------------------------------
E_range = [-1 2];
%---------------------- Load POSCAR ---------------------------------------
[Poscar] = import_poscar('POSCAR');   
coords   = Poscar.coords;
%--------------------------------------------------------------------------

%-------------- Primitive lattice vectors ---------------------------------
PLV  = Poscar.lattice;
%--------------------------------------------------------------------------

%-------------- Primitive reciprocal lattice vectors ----------------------
rec =  inv(PLV)'; 
%--------------------------------------------------------------------------

%% ----------------- Initial info. ----------------------------------------
norb = 4;
Ef   =  0.2;
%% Kpoints %%%
HSP=[  
     0.0000000000   0.0000000000    0.0000000000  %  \Gamma    
     0.0000000000   0.5000000000    0.0000000000  %   X    
     0.5000000000   0.5000000000    0.0000000000  %   M  
     0.0000000000   0.0000000000    0.0000000000  %  \Gamma    
     0.5000000000   0.5000000000    0.5000000000  %   R    
     0.0000000000   0.5000000000    0.0000000000  %   X
 ];



 nk =1000;
[kpoints,nks,list] = KMesh(HSP,nk);

for ikk=1:nks
        
    
        % Transfer the reciprocal to cartesian coordinates.
        kx(ikk)=kpoints(ikk,1)*rec(1,1)+kpoints(ikk,2)*rec(2,1)+kpoints(ikk,3)*rec(3,1);
        ky(ikk)=kpoints(ikk,1)*rec(1,2)+kpoints(ikk,2)*rec(2,2)+kpoints(ikk,3)*rec(3,2);
        kz(ikk)=kpoints(ikk,1)*rec(1,3)+kpoints(ikk,2)*rec(2,3)+kpoints(ikk,3)*rec(3,3);
        % Distance between two kpoints
        if ikk==1           
            kdr(ikk)=0;
            kd(ikk)=0;
        end
        if ikk~=1       
            kdr(ikk)=sqrt((kx(ikk)-kx(ikk-1))^2+(ky(ikk)-ky(ikk-1))^2+(kz(ikk)-kz(ikk-1))^2);  % Distance between two kpoints          
            kd(ikk)=kd(ikk-1)+kdr(ikk);   % add the distance        
        end 
        
end

%% ----------- Eigenvalue and Eigenvector calculations --------------------
tic
Ek     = zeros(nks,norb);

for ik=1:nks
    time_init=tic;
    
    H_SOC = RhSi_WoSOC.RhSi_H(kpoints(ik,:));
    [vec, Etemp] = eig(H_SOC);
    Ek(ik,:)  = diag(Etemp);
    
    if mod(ik,10)==0
        fprintf('%3i/%i: %.3fs [%.2fm]\n',ik,nks,toc,toc(time_init)/60);
    end
    
end
toc

%% ---------------------- Plotting ----------------------------------------
 h = figure('position',[150 0 1200 800]);
hs= size(HSP,1);
kd=kd/2/pi;
hk=nks/(hs);
    for ii=1:norb
        plot(kd,Ek(1:nks,ii)-Ef,'Color',[0.600000023841858 0.200000002980232 0],'LineWidth',2);
        hold on
    end
    
   
%  p1 = plot(kd,Ek(1:nks,29)-Ef,'Color',[0.392156862745098 0.831372549019608 0.0745098039215686],'LineWidth',1.8);
% p2 = plot(kd,Ek(1:nks,2)-Ef,'k-','LineWidth',1.5);
% p3 = plot(kd,Ek(1:nks,3)-Ef,'c-','LineWidth',1.5);
     
 x=[kd(1),kd(end)];
 y=[0 0];
 plot(x,y,'--k')  
        
for hh=1:hs
    x=[kd(list(hh)) kd(list(hh))];
    y=[E_range(1) E_range(2)];
    plot(x,y,'-k')
end
ylabel('\bf{Energy (eV)}','FontSize',24,'interpreter','LaTex');

 axis([kd(1) kd(nks) E_range(1) E_range(2)]);
ax = gca;
ax.FontSize   = 20;
ax.FontWeight = 'bold';
ax.TickLength = [0.02 0.02];
 ax.XTick      = kd(list(:));
 ax.XTickLabel =  {'$\Gamma$','X','M','$\Gamma$','R','X'};
ax.LineWidth  = 2;
ax.TickLabelInterpreter='latex';

 %leg=legend([p1 p2 p3],{'Band-96','Band-97','Band-98'},'interpreter','LaTex');
% leg=legend(p1,{'Nodal-Line Band -- U=0'},'interpreter','LaTex','FontSize',16);
 
 
 
F005  = fullfile(['Band_NSOC'  '.fig']); 
F006  = fullfile(['Band_NSOC'  '.png']);
saveas(h,F005)
saveas(h,F006)
 save Band_Data_YHV.mat 