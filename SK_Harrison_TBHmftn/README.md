# 🧭 SK_Harrison_TBHmftn

This program constructs a ftn58sparse through Slater-Koster-Harrison method.

## [SK method](http://yclept.ucdavis.edu/course/240A.F17/supple/TightBindingPrimer.pdf)

## [Harrison method](./Harrison_ref)

## 🔰 Getting Strated

0️⃣ If you want to use the Harrison method to construct SK model, please uncomment the following line in 'SK_Harrison_Build.m'.

```Matlab
t = slaterkoster(vec2(1,ii)/r2(ii),vec2(2,ii)/r2(ii),vec2(3,ii)/r2(ii),TBK1,TBK2);                                                                         
%t = harrison(TBK1,TBK2,t,r2(ii),powd_1,powd_2,rd(iorb1),rd(iorb2),ro_1,ro_2,ro_3,ro_4);
```

1️⃣ Prepare the structure input file (st_input.txt) for your Wannier tight-binding model. The information should follow the data in .win file.

```txt
SnTe_bulk                          % Name for the model
6.38865 6.38865 6.38865            % Lattice constant  
0.0   0.5   0.5                    % Bravais lattice vectors, where a = (norm(a)/norm(BR(1,:))*(BR(1,1)x + BR(1,2)y + BR(1,3)z) ... etc.
0.5   0.0   0.5
0.5   0.5   0.0
2                                  % Number of element species
Sn 2 3 4                           % (1st type) with the orbit index
Te 2 3 4                           % (2nd type) the order should be the same as part between "begin atoms_cart" and "end atoms_cart" in .win file.
1 1                                % Number of atom for each species, for exampel, 1 Sn and 1 Te. 
 0.00000   0.00000   0.00000    3  % FARCTIONAL COORDINATES for each atom followed by the number of orbits. 
-0.50000   0.50000   0.50000    3  % the order should be the same as part between "begin atoms_cart" and "end atoms_cart" in .win file.
```

,in which the orbit index should follow the table in TBHmftn.m:

```Matlab
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
```

- On the other hand, one can also change a line in TBHmftn.m to match the orbtial labels:

```Matlab
id2txt = {'s' 'px' 'py' 'pz' 'dz' 'dxz' 'dyz' 'dxy' 'dx2y2' ...
          'fz3' 'fxz2' 'fyz2' 'fz(x2-y2)' 'fxyz' 'fx(x2-3y2)' 'fy(3x2-y2)'};
```

2️⃣ Prepare the input in 'SK_Harrison_Build.m'.

``` Matlab
%% --- Input Arguments --- %%%
%%% ---------------------------------------- %%%
fid   = fopen('st_input.txt','r');
Ecut  = 100;                        %--- Energy cutoff   
Rcut  = 5.692;                      %--- Distance cutoff
isSO  = 1;                          %--- switch of Spin-Orbital-coupling  
```

3️⃣ Prepare the SK parameters in 'SK_Harrison_Build.m'.

``` Matlab
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
```

3️⃣ Prepare the SOC parameters in 'SK_Harrison_Build.m'.

```Matlab
%% --- SOC parameters --- %%%
lamda_W  =[0.58];  % atomic SOC parameter for d-orbitals
lamda_Se =[0.43];  % atomic SOC parameter for p-orbitals
ipd=[1,6];         % the index of dz2, and for the order dz2,dxz,dyz,dxy,dx2y2
ipp=[11,14,17,20]; % the index of px, and for the order px,py,pz
```

3️⃣ Prepare the onsite energies in 'SK_Harrison_Build.m'.

```Matlab
%% --- Crystal Field (*onsite energy*) --- %%
%onsite        = [1:norb];
onsite(1:5)   = [3.93087;5.45;5.45;4.50;4.50]; 
onsite(6:10)  = [3.93087;5.45;5.45;4.50;4.50]; 
onsite(11:13) = [3.43;3.43;2.70741]; 
onsite(14:16) = [3.43;3.43;2.70741];
onsite(17:19) = [3.43;3.43;2.70741];
onsite(20:22) = [3.43;3.43;2.70741];   
```

3️⃣ Prepare the parameters for Harrison method in 'SK_Harrison_Build.m' if you want to use Harrison method to construct SK model.

```Matlab
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
```

4️⃣ Run 'SK_Harrison_Build.m' to obtain the SKftn58sparse.

```Matlab
>> SK_Harrison_Build.m
```

5️⃣ Run 'bandSK.m' to check the band structure of the SKftn58sparse.

```Matlab
>> bandSK(0,Erange,Ef);
```
