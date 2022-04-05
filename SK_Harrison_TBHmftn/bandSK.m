function bandSK(isSys,E,Ef)
%%% Function for plotting band structure from SKftn58sparse %%%
%%% ------------------------------------------------------- %%%
%%% Input arguments                                         %%%
%%% isSys  => crude SK model or parametrize by ftn58        %%%
%%%           (0 for crude SK; 1 for ftn58)                 %%%
%%% ------------------------------------------------------- %%%

Eup = E(2);
Edn = E(1);

%% --- Load SKftn58sparse data and kpoints --- %%%
if isSys==1
    load DFTSKftn58sparse.mat
    norb = DFTSKftn58sparse.norb;
    ii   = DFTSKftn58sparse.ij(:,1);
    jj   = DFTSKftn58sparse.ij(:,2);
    dd   = DFTSKftn58sparse.dd;
    tt   = DFTSKftn58sparse.tt; 
else
    load SKftn58sparse.mat
    norb = SKftn58sparse.norb;
    ii   = SKftn58sparse.ij(:,1);
    jj   = SKftn58sparse.ij(:,2);
    dd   = SKftn58sparse.dd;
    tt   = SKftn58sparse.tt;
end
    
    nk  = 50;
    p1 = [linspace(0.0,0.5,nk+1)' linspace(0.0,0.0,nk+1)' linspace(0.0,0.0,nk+1)'];
    p2 = [linspace(0.5,1/3,nk+1)' linspace(0.0,1/3,nk+1)' linspace(0.0,0.0,nk+1)'];
    p3 = [linspace(1/3,0.0,nk+1)' linspace(1/3,0.0,nk+1)' linspace(0.0,0.0,nk+1)'];
    p4 = [linspace(0.0,0.0,nk)' linspace(0.0,0.0,nk)' linspace(0.0,0.5,nk)'];
    
    kpoints = [p1(1:nk,:);p2(1:nk,:);p3(1:nk,:);p4(1:nk,:)];

%% --- Plot band structrue --- %%%
tic
eigvec = cell(size(kpoints,1),1);
for ik=1:size(kpoints,1)
    kcolumnvec=kpoints(ik,:)'*2*pi;
    Hsparse=sparse(ii,jj,exp(1i*dd*kcolumnvec).*tt,norb,norb);
    HH=full(Hsparse);
    HH=(HH+HH')/2;
    [eigvectemp ,Etemp]=eig(HH);
     
    eigvec{ik,1} = eigvectemp;

    A(ik*norb-norb+1:ik*norb,1)=ik;
    A(ik*norb-norb+1:ik*norb,2)=diag(Etemp);
    
    Ek_sk(ik,:) = diag(Etemp);    
end
toc

%% --- Plotting Details --- %%%
figure('position',[150 0 850 660],'paperposition',[0.25 0.25 8 10.5],...
       'papertype','usletter','numbertitle','off',...
       'PaperPositionMode','manual','paperorientation','portrait',...
       'color','w');

nks = length(kpoints);
% for ii=1:fnorb
%     h1 = plot([1:nks],Ek(:,ii)-Ef,'b-','MarkerSize',5);
%     hold on
% end

for ii=1:norb
    h1 = plot([1:nks],Ek_sk(:,ii)-Ef,'r-','MarkerSize',5);
    hold on
end

axis([A(1,1) A(end,1) Edn Eup]);
ylabel('\bf{Energy (eV)}','interpreter','LaTex');
ls = legend(h1,'\bf{SK}','Location','best');
set(ls,'interpreter','LaTex','FontSize',18);
line('XData', [0 A(end,1)], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 0.5, 'Color','k');

list=[1 51 101 151 200];
for il = 1:size(list,2)
line('XData', [list(il) list(il)], 'YData', [Edn Eup], 'LineStyle', '-', ...
    'LineWidth', 0.1, 'Color','k');
end

ax = gca;
ax.FontSize = 24;
ax.FontWeight = 'bold';
ax.TickLength = [0.02 0.02];
ax.XTick = list(:);
ax.TickLabelInterpreter='latex';
ax.XTickLabel = {'$\Gamma$' 'M' 'K' '$\Gamma$','Z'};
ax.LineWidth = 1.2;

save eigvec.mat eigvec