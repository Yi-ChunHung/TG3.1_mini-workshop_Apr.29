clear all

%%% --- inputs --- %%%
load Narrow_Ribbon_ftn58sparse.mat
filling    = 2;       % number of filled bands on each site
nocc       = 3200;
kvec       = 0; % at kvec
minbondlen = 2.6;     % only bonds with bond length > minbondlen would be plot
maxbondlen = 3;       % only bonds with bond length < maxbondlen would be plot
np         = 16;

SPftn58sparse = NRftn58sparse;
%%% --- Initialization --- %%%
ii    = SPftn58sparse.ij(:,1);
jj    = SPftn58sparse.ij(:,2);
tt    = SPftn58sparse.tt;
dd    = SPftn58sparse.dd;
natom = SPftn58sparse.Nat;
norb  = SPftn58sparse.norb;
BR    = SPftn58sparse.BR1D;
abc   = SPftn58sparse.abc;
BR    = diag(abc)*BR;
clearvars abc NRftn58sparse

%%% --- Real Procedure --- %%%
figure
set(gcf, 'unit', 'normalized', 'position', [0.4, 0.5, 0.38, 0.68])
set(gca, 'Outerposition', [0.0, 0.0, 1, 1], 'Position', [0.01, 0.03, 0.98, 0.94])

%% Plotting bonds
fprintf('Start plotting bonds ... \n')
tic
apos = zeros(natom,3);
for i = 1:natom
    apos(i,:) = SPftn58sparse.Ainfo(i).Position*BR;
end
bondtable     = nchoosek(1:natom,2);
for i = 1:length(bondtable)
    if (norm(apos(bondtable(i,1),:) - apos(bondtable(i,2),:)) > minbondlen  & norm(apos(bondtable(i,1),:) - apos(bondtable(i,2),:)) < maxbondlen)
        line([apos(bondtable(i,1),1) apos(bondtable(i,2),1)], [apos(bondtable(i,1),2) apos(bondtable(i,2),2)], [apos(bondtable(i,1),3) apos(bondtable(i,2),3)], 'linew', 1.5,...
            'color', [0.8 0.8 0.8]); % [0.6 0.75 0.75]
    end
end
clearvars BR bondtable apos
toc

%% Calculate the eigenvectors
fprintf('Start calculating eigens ... \n')
tic
HH        = full(sparse(ii,jj,tt.*exp(1i*dd*2*pi*kvec'),norb,norb));
HH        = (HH + HH')/2;
[efun, E] = eig(HH);
E         = diag(E);
[~,occ]   = sort(E);
save Eigen.mat efun E -v7.3 
clearvars E HH
toc

V = efun(:,occ(1:nocc));

%% Construct table of nano-disk sites
Orbitps       = SPftn58sparse.Orbitps;
orbpos        = Orbitps(:,4:6);
[posxy,~,icc] = unique(orbpos(:,1:2),'rows'); % orbpos = posxy(icc,:)
[uicc,~,~]    = unique(icc);                  
Orbitps       = [Orbitps,icc];
delta_rho     = zeros(size(posxy,1),1);
delta_rho_t   = zeros(size(posxy,1),1);
clearvars icc orbpos

%% Calculate the charge distribution
c = parcluster('local');
c.NumWorkers = np;
parpool(c, c.NumWorkers);

tic
fprintf('Start calculating PCD ... \n')
parfor ic = 1:length(uicc)
    orbid           = find(Orbitps(:,end)==uicc(ic));
    delta_rho_t(ic) = sum(abs(V(orbid,:)).^2,'all') - filling;
end
delta_rho = delta_rho_t(uicc);

clearvars proj orbid delta_rho_t uicc V
save delta_rho.mat delta_rho posxy -v7.3
toc

%% Plot the RSPCD
fprintf('Start plotting PCD ... \n')
tic
hold on
scatter(posxy(:, 1), posxy(:, 2), 10, delta_rho, 'filled')
colorbar
toc
axis off
axis equal

sum(delta_rho)

delete(gcp)