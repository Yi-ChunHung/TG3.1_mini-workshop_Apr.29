%% Author: Bao-Kai, Wang
%% Editor: Yi-Chun, Hung
clear all

%%% --- inputs --- %%%
load Narrow_Ribbon_ftn58sparse.mat
state      = 109;     % the state^th eigenvector
kvec       = [0 0 0]; % at kvec
minbondlen = 2.6;     % only bonds with bond length > minbondlen would be plot
maxbondlen = 3;       % only bonds with bond length < maxbondlen would be plot

SPftn58sparse = NRftn58sparse;
%%% --- Initialization --- %%%
ii    = SPftn58sparse.ij(:,1);
jj    = SPftn58sparse.ij(:,2);
tt    = SPftn58sparse.tt;
dd    = SPftn58sparse.dd;
natom = SPftn58sparse.Nat;
norb  = SPftn58sparse.norb/2;
BR    = SPftn58sparse.BR1D;
abc   = SPftn58sparse.abc;
BR    = diag(abc)*BR;

%%% --- Real Procedure --- %%%
figure
set(gcf, 'unit', 'normalized', 'position', [0.4, 0.5, 0.38, 0.68])
set(gca, 'Outerposition', [0.0, 0.0, 1, 1], 'Position', [0.01, 0.03, 0.98, 0.94])

%% Plot the bonds
fprintf('Start plotting bonds ... \n')
tic
%for ii = 1 : natom - 1
%    for jj = ii+1 : natom
%        apos_ii  = SPftn58sparse.Ainfo(ii).Position*BR;
%        apos_jj  = SPftn58sparse.Ainfo(jj).Position*BR;
%        if (norm(apos_ii - apos_jj) > minbondlen  & norm(apos_ii - apos_jj) < maxbondlen)
%            line([apos_ii(1) apos_jj(1)], [apos_ii(2) apos_jj(2)], [apos_ii(3) apos_jj(3)], 'linew', 1.5,...
%                'color', [0.8 0.8 0.8]); % [0.6 0.75 0.75]
%        end
%    end
%end
apos = zeros(natom,3);
for i = 1:natom
    apos(i,:) = SPftn58sparse.Ainfo(i).Position*BR;
end
bondtable = nchoosek(1:natom,2);
for i = 1:length(bondtable)
    if (norm(apos(bondtable(i,1),:) - apos(bondtable(i,2),:)) > minbondlen  & norm(apos(bondtable(i,1),:) - apos(bondtable(i,2),:)) < maxbondlen)
        line([apos(bondtable(i,1),1) apos(bondtable(i,2),1)], [apos(bondtable(i,1),2) apos(bondtable(i,2),2)], [apos(bondtable(i,1),3) apos(bondtable(i,2),3)], 'linew', 1.5,...
            'color', [0.8 0.8 0.8]); % [0.6 0.75 0.75]
    end
end
toc

%% Calculate the eigenvectors
fprintf('Start calculating eigens ... \n')
tic
HH        = full(sparse(ii,jj,tt.*exp(1i*dd*2*pi*kvec'),norb,norb));
HH        = (HH + HH')/2;
[efun, ~] = eig(HH);
toc

%% Calculate the RSPCD
fprintf('Start calculating PCD ... \n')
tic
for iatom = 1:natom
    orbid         = SPftn58sparse.Ainfo(iatom).OrbitID;
    orbspec       = length(orbid);  % how many (spin up) orbitals on an atom   
    spinup(iatom) = sum(abs(efun((iatom-1)*orbspec+1 : iatom*orbspec, state)).^2);
    spindn(iatom) = sum(abs(efun((iatom-1)*orbspec+1+norb: iatom*orbspec+norb, state)).^2);
    %spinup(iatom) = efun((iatom-1)*orbspec+1 : iatom*orbspec, state)'*eye(orbspec)*efun((iatom-1)*orbspec+1 : iatom*orbspec, state);
    %spindn(iatom) = efun((iatom-1)*orbspec+1+norb : iatom*orbspec+norb, state)'*eye(orbspec)*efun((iatom-1)*orbspec+1+norb : iatom*orbspec+norb, state);
end
toc

%% Plot the RSPCD
fprintf('Start plotting PCD ... \n')
tic
hold on
scatter3(apos(:, 1), apos(:, 2), apos(:, 3), 15000*(spinup(:) + spindn(:)), [0.0 0.8 0], 'filled')
toc
axis off
axis equal
