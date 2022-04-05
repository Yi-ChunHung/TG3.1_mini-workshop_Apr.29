%%%   Program for buildig ftn58 based on the Slater-Koster model     %%%
%%% ---------------------------------------------------------------- %%%
%%% Necessary Input files:                                           %%%
%%% 1) st_input.txt                                                  %%%
%%% ---------------------------------------------------------------- %%%
%%% Required Functions:                                              %%%
%%% 1) latticegen.m                                                  %%%
%%% 2) slaterkoster.m                                                %%%
%%% ---------------------------------------------------------------- %%%

%%% ---------------------------------------------------------------------- %%%
%%%                     Order of atomic orbitals: 
%%%    id   Wien2k       QE(w90)      OpenMX       TB2KKR
%%% s   1   s            s            s            s
%%% p   2   px           pz           px           px
%%%     3   py           px           py           py
%%%     4   pz           py           pz           pz
%%% d   5   dxy          dz2          dz2          dz2
%%%     6   dyz          dxz          d(x2-y2)     dxz
%%%     7   dxz          dyz          dxy          dyz
%%%     8   d(x2-y2)     d(x2-y2)     dxz          dxy
%%%     9   dz2          dxy          dyz          d(x2-y2)
%%% f  10   fy(3x2-y2)   fz3          fz3          fz3        (z(5z2-3r2))
%%%    11   fz(x2-y2)    fxz2         fxz2         fxz2       (x(5z2-3r2))  
%%%    12   fyz2         fyz2         fyz2         fyz2       (y(5z2-3r2))
%%%    13   fxz2         fz(x2-y2)    fz(x2-y2)    fz(x2-y2)  (z(x2-y2))
%%%    14   fxyz         fxyz         fxyz         fxyz       (2xyz)
%%%    15   fx(x2-3y2)   fx(x2-3y2)   fx(x2-3y2)   fx(x2-3y2) (x(x2-3y2))
%%%    16   fz3          fy(3x2-y2)   fy(3x2-y2)   fy(3x2-y2) (y(3x2-y2))
%%% --------------------------------------------------------------------- %%%
clear all

%% --- Input Arguments --- %%%
%%% ---------------------------------------- %%%
fid   = fopen('st_input.txt','r');
Ecut  = 100;                        %--- Energy cutoff   
Rcut  = 5.692;                      %--- Distance cutoff
isSO  = 1;                          %--- switch of Spin-Orbital-coupling   

%% --- SK Parameters (first neighbor) --- %%%
global sss sps pps ppp sds pds pdp dds ddp ddd
sss = -0.71187;                                    
sps =  0.70334;
pps =  2.35019;
ppp = -1.006;
sds = -3.12;
pds = -4.82;
pdp =  4.08;
dds = -4.22;
ddp =  3.60;
ddd = -2.29;
%%% ---------- %%%

%% --- SOC parameters --- %%%
lamda_W  =[0.58];%0.58
lamda_Se =[0.43];%0.43
ipd=[1,6];
ipp=[11,14,17,20];

%% --- Crystal Field (*onsite energy*) --- %%
%onsite        = [1:norb];
onsite(1:5)   = [3.93087;5.45;5.45;4.50;4.50]; 
onsite(6:10)  = [3.93087;5.45;5.45;4.50;4.50]; 
onsite(11:13) = [3.43;3.43;2.70741]; 
onsite(14:16) = [3.43;3.43;2.70741];
onsite(17:19) = [3.43;3.43;2.70741];
onsite(20:22) = [3.43;3.43;2.70741];                

%% --- Adjusted distance dependent TB parameters (Harrison) --- %%
%powd_1 = 4.4;       % order of atoms' distance for s,p & d
%powd_2 = 4.4;       % order of atoms' distance for d & d
%ro_1 = 1;           % order of rd(orb1) , s,p & d
%ro_2 = 1;           % order of rd(orb2) , s,p & d
%ro_3 = 1;           % order of rd(orb1) , d & d
%ro_4 = 1;           % order of rd(orb2) , d & d

%rd = ones([1 norb]);
%rd([1 6]) = 1.0;
%rd([2:3 7:8]) = 1.0;
%rd([4:5 9:10]) = 1.0;
%rd([11:12 14:15 17:18 20:21]) = 1.0;
%rd([13 16 19 22]) = 1.0;
%%% ---------------------------------------- %%%
%%% ---------------------------------------- %%%

%%% --------------------------------------------------------------------- %%%
%% --- Reading information from st_input.txt to construct ftn58sparse --- %%%
buf    = fgetl(fid);                    % system's name
system = buf;
buf    = fgetl(fid);                    % lattice const.
ABC    = str2num(buf);

for ii=1:3                              % get position info.
    buf = fgetl(fid);
    BR_1(ii,:) = str2num(buf);          % BR_1(ii,:) ii:1->3 
end

buf = fgetl(fid);                       % fifth line --> the sustem's atom type
atn = str2num(buf);                     % # of atom type

orbit = zeros(atn,5);
for ii=1:atn                            
    buf        = fgetl(fid);        
    space      = find(buf(:)==' ');                                         
                                                                            
        %fprintf('%d\n',space);
    norb       = length(space);                                             
        %fprintf('%d\n',norb);
    atom(ii,:) = buf(1:(space(1)-1));                                       % atom name
    orbit(ii,1:length(space)) = str2num(buf((space(1)+1):end));             % which orbitals (id)
end                                     

buf = fgetl(fid);
nat = str2num(buf);                     % # of atom
%fprintf('%d\n',length(nat));

id2txt = {'s' 'px' 'py' 'pz' 'dz' 'dxz' 'dyz' 'dxy' 'dx2y2' 'f-3' 'f-2' 'f-1' 'f0' 'f1' 'f2' 'f3'};

id_at   = 1;                            
id_orb  = 1;                            
id_orbp = 1;                            
for aa=1:length(nat)                    % Type of atoms
    for bb=1:nat(aa)                    % # of atoms
        buf                             = fgetl(fid);
        ps                              = str2num(buf);
        atominfo(id_at).Atom(1,:)       = atom(aa,:);                       % atom name
        atominfo(id_at).Position(1,1:3) = ps(1:3);                          % atom position
        atominfo(id_at).Norb            = ps(4);                            % # of atom's orbital
        atominfo(id_at).OrbitIndex      = id_orb:(id_orb+ps(4)-1);          % give orbital index
        
        temob    = [];
        orbit_n0 = find(orbit(aa,:)>0);                                     
        for kk=1:length(orbit_n0)                                           % label which orbital
            temob = strcat(temob,{' '},id2txt(1,orbit(aa,kk)));
        end
        atominfo(id_at).Orbit(1,:) = cellstr(temob);
        
        %%% Orbitps format:
        %%% #_of_orbital orbital_index(ftn58) which_atom frac_coord standard_orbital_id  
        for kk=1:ps(4)
            orbitps(id_orbp,1:7) = [id_orbp,id_orb+(kk-1),id_at,ps(1:3),orbit(aa,kk)];
            id_orbp              = id_orbp + 1;
        end
        
        atominfo(id_at).OrbitID  = orbit(aa,1:length(orbit_n0));
        id_at                    = id_at + 1;
        id_orb                   = id_orb + ps(4);
    end 
end
SKftn58sparse.System  = system;
SKftn58sparse.Ainfo   = atominfo;
SKftn58sparse.abc     = ABC;
SKftn58sparse.BR      = BR_1;               
SKftn58sparse.Nat     = length(atominfo);
SKftn58sparse.ver     = 'type1';
SKftn58sparse.norb    = size(orbitps,1);
SKftn58sparse.Orbitps = orbitps;

%%% ---------------------------------------------- %%%
%%% ---------------------------------------------- %%%
%% --- Actual Procedure for Building up SKftn58 --- %%%
abc     = SKftn58sparse.abc;
BR      = SKftn58sparse.BR;
norb    = SKftn58sparse.norb;                   
%abc2xyz = diag(abc)*BR;

%--- real space position of cells (hopping distance)
lat_BR   = latticegen(BR,5,5,5);            % cell position in x,y,z basis
lat_unit = latticegen(eye(3),5,5,5);        % cell position in a,b,c basis              
nlat     = size(lat_BR,2);                  % size(A,2) : # of column       

%--- TB2KKR Format 
%--- orbit_id orbit_atom_id global_orbit_id(see the Table)
TB2KKR = [orbitps(:,1) orbitps(:,3) orbitps(:,7) ];

%% --- Orbital positions in Cartesian basis --- %%
pos = (orbitps(:,4:6)*BR)';

tic
%ftn58 = [];
ftn58 = cell(norb*(norb+1)/2,1);
for iorb1=1:norb
    pos1 = pos(:,iorb1);
    for iorb2 = iorb1:norb
        
        %--- pos2: positions of the 2nd orbit relative to the 1st orbit  
        pos2   = repmat(pos(:,iorb2)-pos1,1,nlat)+lat_BR;
        vec2   = diag(abc)*pos2;  
        %--- r2: distances of each intended hopping path
        r2     = sqrt(sum(vec2.*vec2,1));      
        iself  = find(r2<1e-4);         %--- hopping within a same atom
        iself  = [iself find(r2>Rcut)]; %--- Look for unnecessary hopping atoms 
        hop_id = setdiff(1:nlat,iself);
         
        %ftn58temp = [];
        ftn58temp = cell(length(hop_id)+1,1);
        for ii = hop_id
            %HOP(iorb1,iorb2,ii) = ii;
            %--- hoppings from irob1 to irob2 
            TBK1 = TB2KKR(iorb1,3);
            TBK2 = TB2KKR(iorb2,3);
            
            t = slaterkoster(vec2(1,ii)/r2(ii),vec2(2,ii)/r2(ii),vec2(3,ii)/r2(ii),TBK1,TBK2);                                                                         
            %t = harrison(TBK1,TBK2,t,r2(ii),powd_1,powd_2,rd(iorb1),rd(iorb2),ro_1,ro_2,ro_3,ro_4);
            
            if abs(t) <= Ecut
                %ftn58temp = [ftn58temp;[iorb1 iorb2 t lat_unit(:,ii)']];
                ftn58temp{ii} = [iorb1 iorb2 t lat_unit(:,ii)'];
            end
        end
        if iorb1==iorb2
            %ftn58temp = [ftn58temp; [iorb1 iorb2 onsite(iorb1) 0 0 0]];
            ftn58temp{length(hop_id)+1} = [iorb1 iorb2 onsite(iorb1) 0 0 0];
        end
        %ftn58 = [ftn58;ftn58temp];
        ftn58{(iorb1-1)*norb+((iorb2-iorb1)+1)} = cell2mat(ftn58temp);
    end
end
ftn58 = cell2mat(ftn58);
toc

nbond  = size(ftn58,1);
ftn58  = [(1:nbond)' ftn58];
ftn58  = [[norb nbond 0 0 0 0 0];ftn58];
nn_hop = max(ftn58(2:end,5:7));
fprintf('Max Hopping Site = %3i\n',nn_hop);

%save TBK.mat TBK TBK_typ HOP;

%% --- Save Info. to SKftn58sparse.mat --- %%%
if isSO==0
    SKftn58sparse.isSO = 0;
    ib1=find(ftn58(2:end,2) == ftn58(2:end,3))+1;
    ib2=find(ftn58(2:end,2) ~= ftn58(2:end,3))+1;
    SKftn58sparse.norb       = ftn58(1,1);
    SKftn58sparse.Orbitps    = orbitps;
    SKftn58sparse.ij         = [ftn58(ib1,2:3);ftn58(ib2,2:3);ftn58(ib2,[3 2])];
    SKftn58sparse.tt         = [ftn58(ib1,4);ftn58(ib2,4);ftn58(ib2,4)];
    SKftn58sparse.dd         = [ftn58(ib1,5:7);ftn58(ib2,5:7);-ftn58(ib2,5:7)];
else
    [ftn58, ftn58SO] = SOftn58(0,ftn58,lamda_Se,ipp,lamda_W,ipd);

    SKftn58sparse.isSO = 1;
    ib1=find(ftn58(2:end,2) == ftn58(2:end,3))+1;
    ib2=find(ftn58(2:end,2) ~= ftn58(2:end,3))+1;
    SKftn58sparse.norb       = ftn58(1,1);
    SKftn58sparse.Orbitps    = [orbitps;orbitps];
    SKftn58sparse.ij         = [ftn58SO(2:end,2:3);ftn58(ib1,2:3);ftn58(ib2,2:3);ftn58SO(2:end,[3 2]);ftn58(ib2,[3 2])];
    SKftn58sparse.tt         = [ftn58SO(2:end,4);ftn58(ib1,4);ftn58(ib2,4);conj(ftn58SO(2:end,4));ftn58(ib2,4)];
    SKftn58sparse.dd         = [ftn58SO(2:end,5:7);ftn58(ib1,5:7);ftn58(ib2,5:7);-ftn58SO(2:end,5:7);-ftn58(ib2,5:7)];    
end
save SKftn58sparse.mat SKftn58sparse

%%% Check the band structure of the SK model
%Erange = [-15 5];
%Ef     = 6.5;
%bandSK(0,Erange,Ef);
%toc

%whos SKftn58sparse
