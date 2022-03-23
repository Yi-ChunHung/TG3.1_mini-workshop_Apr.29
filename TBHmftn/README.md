# 🧭 TBHmftn

This script will convert the case_hr.dat generated from the wannier90 or a ftn58 to a packed file named ftn58sparse.mat.

## 🔰 Getting Strated

0️⃣ Choose the source for ftn58sparse.mat in TBHmftn.m

```Matlab
%% --- Choose your ftn58 resource --- %%
if 1==1
    [ftn58,~]=readwanhr(wcal.ref,0);   
%     save ftn58.mat
else 
    load ftn58.mat
end
```

1️⃣ Prepare the structure input file (st_input.txt) for your Wannier tight-binding model. The information should follow the data in .win file. An example for SnTe is given below:

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

2️⃣ Prepare the input (input.txt) file for the program "TBHmftn".

``` txt
TBHmftn                         

wcal.isEC   = 0                % isEC = 1 to truncate the TB model based on the size of hopping strength
wcal.ecut   = 0                % specify the truncated hopping strength
wcal.isSO   = 1                % isSO = 0 model w/o SOC; isSO = 1 model w/ SOC
wcal.ref    = SnTeso_hr.dat    % name for the TB model generated from Wannier90.x 

endTBHmftn
```

3️⃣ Run "TBHmftn" to get the packed file "ftn58sparse.mat"

```Matlab
>> TBHmftn.m
```

4️⃣ Get the benchmark for your wannier TB model.

```Matlab
>> band_Bench.m
```

5️⃣ Calculate the band structure along the desired k paths

```Matlab
>> band_ftn.m
```
