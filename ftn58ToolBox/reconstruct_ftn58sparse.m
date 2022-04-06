function rc_ftn58sparse = reconstruct_ftn58sparse(SPftn58sparse,removed_atoms,name)
% The program remove certain atoms in ftn58sparse. Then, it reconstruct ftn58sparse based on this change
% SPftn58sparse is not necessarily to be the supercell ftn58saprse
% removed_atoms = [5,6,10,12]; ... etc.

rc_ftn58sparse        = SPftn58sparse;
rc_ftn58sparse.System = name;
isSO                  = rc_ftn58sparse.isSO;
norb_origin           = SPftn58sparse.norb;   
rc_ftn58sparse.Nat    = rc_ftn58sparse.Nat - length(removed_atoms);

%% reconstructed Orbitps
% ----------------------------------------------------------------------------------
% Format of orbitps 
% 1     => orbital ID 
% 2,3   => original orbital ID and atom ID in ftn58sparse
% 4:6   => orbital frac. coord. 
% 10    => orbital standard ID (see table)
% SPftn58sparse.Orbitps = [orbitps(1:end,1:3) orbitps(1:end,4:6) orbitps(1:end,10)];
% ----------------------------------------------------------------------------------
Orbitps          = rc_ftn58sparse.Orbitps;
removed_orbitals = [];                                         % codes for reconstructed ij, dd, tt
% column 3
for i = 1:length(removed_atoms)
    atom = removed_atoms(i);
    Orbitps(find(Orbitps(:,3)==atom),:) = [];

    orbitals           = rc_ftn58sparse.Ainfo(atom).OrbitIndex; % codes for reconstructed ij, dd, tt
    removed_orbitals   = [removed_orbitals,orbitals];           % codes for reconstructed ij, dd, tt
end
if isSO
    removed_orbitals   = [removed_orbitals,removed_orbitals+(norb_origin/2)]; % codes for reconstructed ij, dd, tt
end
atomID_temp            = unique(Orbitps(:,3));
[~,atomID]             = ismember(Orbitps(:,3),atomID_temp);
Orbitps(:,3)           = atomID;
% column 1
norb                   = size(Orbitps,1);  
rc_ftn58sparse.norb    = norb;
Origin_OrbitID         = Orbitps(:,1);   % domain of the map from original OrbitalID to the reindexed one
Orbitps(:,1)           = 1:norb;         % reindex the OrbitalID
rc_ftn58sparse.Orbitps = Orbitps;

%% reconstructed ij, dd, tt
indices = find( ismember(rc_ftn58sparse.ij(:,1),removed_orbitals) | ismember(rc_ftn58sparse.ij(:,2),removed_orbitals) );
rc_ftn58sparse.ij(indices,:) = [];
rc_ftn58sparse.tt(indices)   = [];
rc_ftn58sparse.dd(indices,:) = [];
% map the original OrbitalID to the reindexed one
[~,ii]            = ismember(rc_ftn58sparse.ij(:,1),Origin_OrbitID);
[~,jj]            = ismember(rc_ftn58sparse.ij(:,2),Origin_OrbitID);
rc_ftn58sparse.ij = [ii,jj];

%% reconstructed Ainfo
rc_ftn58sparse.Ainfo(removed_atoms) = [];
for atom = 1:rc_ftn58sparse.Nat
    OrbitIndex = Orbitps(find(Orbitps(:,3)==atom),1)';
    if isSO
        OrbitIndex(OrbitIndex>(norb/2)) = [];           % neglect the other spin just as the convention
    end
    rc_ftn58sparse.Ainfo(atom).OrbitIndex = OrbitIndex;
end
