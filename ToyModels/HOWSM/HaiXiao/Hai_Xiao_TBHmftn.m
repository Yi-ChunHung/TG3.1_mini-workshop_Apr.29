%%% Hai-Xiao's Model
%%% arXiv:2007.05068v1 [cond-mat.mtrl-sci] 9 Jul 2020

%%--- Parameters --- %%
onsite = 0;
% NN hoppings
gamma  = 0.6;
lambda = 1;
% NNN hoppings
t1     = 0.2;
t2     = 0.24;

%%--- TBmodel --- %%
ftn58sparse.System   = 'Hai-Xiao Model';
ftn58sparse.toymodel = 1;       % quantative description of crystal
ftn58sparse.abc 	 = [1 1 2]; % tetragonal crystal
ftn58sparse.BR  	 = eye(3);
ftn58sparse.isSO     = 0;
ftn58sparse.Nat      = 4;
ftn58sparse.norb     = 4;

%% quantative description to preserve S4z, C2x, C2y symmetries
% origin lies in the center of cell
x = 1/4;
y = 1/4;
z = -3/8; % C2x(1) == 4, C2y(1) == 3
ftn58sparse.Orbitps             = [[1,1,1,x,y,z,1];[2,2,2,-x,-y,z+1/2,1];[3,3,3,y,-x,z+1/4,1];[4,4,4,-y,x,z+3/4,1]];
ftn58sparse.Ainfo(1).Atom       = '1';
ftn58sparse.Ainfo(1).Position   = [x y z];
ftn58sparse.Ainfo(1).Norb       = 1;
ftn58sparse.Ainfo(1).OrbitIndex = 1;
ftn58sparse.Ainfo(1).Orbit      = 's';
ftn58sparse.Ainfo(1).OrbitID    = 1;

ftn58sparse.Ainfo(2).Atom       = '2';
ftn58sparse.Ainfo(2).Position   = [-x -y z+1/2];
ftn58sparse.Ainfo(2).Norb       = 1;
ftn58sparse.Ainfo(2).OrbitIndex = 2;
ftn58sparse.Ainfo(2).Orbit      = 's';
ftn58sparse.Ainfo(2).OrbitID    = 1;

ftn58sparse.Ainfo(3).Atom       = '3';
ftn58sparse.Ainfo(3).Position   = [y -x z+1];
ftn58sparse.Ainfo(3).Norb       = 1;
ftn58sparse.Ainfo(3).OrbitIndex = 3;
ftn58sparse.Ainfo(3).Orbit      = 's';
ftn58sparse.Ainfo(3).OrbitID    = 1;

ftn58sparse.Ainfo(4).Atom       = '4';
ftn58sparse.Ainfo(4).Position   = [-y x z+3/4];
ftn58sparse.Ainfo(4).Norb       = 1;
ftn58sparse.Ainfo(4).OrbitIndex = 4;
ftn58sparse.Ainfo(4).Orbit      = 's';
ftn58sparse.Ainfo(4).OrbitID    = 1;

%ftn58sparse.Ainfo.Atom       = ['1';'2';'3';'4'];
%ftn58sparse.Ainfo.Position   = [[x y z];[-x -y z+1/2];[y -x z+1];[-y x z+3/4]];
%ftn58sparse.Ainfo.Norb       = [1;1;1;1];
%ftn58sparse.Ainfo.OrbitIndex = [1;2;3;4];
%ftn58sparse.Ainfo.Orbit      = ['s';'s';'s';'s'];
%ftn58sparse.Ainfo.OrbitID    = [1;1;1;1];

ij = [];
tt = [];
dd = [];

%% hoppings: eq.(1) & Fig.2 (a) w/ legend in Fig.2 (b)
% onsite
ij  = [ij;[1 1];[2 2];[3 3];[4 4]];
tt  = [tt;onsite;onsite;onsite;onsite];
dd  = [dd;[0 0 0];[0 0 0];[0 0 0];[0 0 0]];

% NN
% intra
ij = [ij;[3 1];[4 2];[3 2]];
tt = [tt;gamma;gamma;gamma];
dd = [dd;[0 0 0];[0 0 0];[0 0 0]];

ij = [ij;[1 3];[2 4];[2 3]];
tt = [tt;gamma;gamma;gamma];
dd = [dd;[0 0 0];[0 0 0];[0 0 0]];

% inter
ij = [ij;[1 4];[3 1];[4 2];[3 2]];
tt = [tt;lambda;lambda;lambda;lambda];
dd = [dd;[0 -1 1];[-1 0 0];[1 0 0];[0 1 0]];

ij = [ij;[4 1];[1 3];[2 4];[2 3]];
tt = [tt;lambda;lambda;lambda;lambda];
dd = [dd;[0 1 -1];[1 0 0];[-1 0 0];[0 -1 0]];

% NNN 
% intra
ij = [ij;[1 2];[3 4]];
tt = [tt; t1; t1];
dd = [dd;[0 0 0];[0 0 0]];

ij = [ij;[2 1];[4 3]];
tt = [tt; t1; t1];
dd = [dd;[0 0 0];[0 0 0]];

% inter
ij = [ij;[1 2];[3 4]];
tt = [tt; t2; t2];
dd = [dd;[1 1 -1];[-1 1 -1]];

ij = [ij;[2 1];[4 3]];
tt = [tt; t2; t2];
dd = [dd;[-1 -1 1];[1 -1 1]];

ftn58sparse.ij 		 = ij;
ftn58sparse.tt       = tt;
ftn58sparse.dd 		 = dd;

save Hai_Xiao_TBHmftn.mat ftn58sparse
