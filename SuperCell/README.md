# 🧭 SuperCell

The program can generate an arbitrary supercell from a TB model with new unit vectors constructed from the linear combination of the original unit vecotrs.

## :notebook: [Principle](./Principle%20of%20supercell.pdf)

## ❗ Warning

The format of Orbitps will greatly affect the validity of this program. Please pay attention on it for the ftn58sparse that may not comes from the standard process.

## 🔰 Getting Started

0️⃣ Check the order of the atoms ('atomps') is sorted in the desired way through the lines at the bottom of 'SuperCell.m'. If there is no special needs for order of the atoms, comment these lines.

```Matlab
%%% --------------------------------------------------- %%%
%%% --- sort the atomps (for Green ftn calculation) --- %%%
%%% --------------------------------------------------- %%%
frac_cord                   = 1:3;
frac_cord(ribbon_direction) = [];
% atomps                     = sortrows(atomps,frac_cord,'descend');
%[~,index]                   = sortrows(floor(atomps),frac_cord,'descend');
%atomps(1:end)               = atomps(index,:);

f1           = frac_cord(1);
f2           = frac_cord(2);
tru_atomps   = atomps;
tru_atomps   = floor(tru_atomps);

NR_11_id     = find(tru_atomps(:,f1)<(nv(f1)/2) & tru_atomps(:,f2)<(nv(f2)/2) );
NR_12_id     = find(tru_atomps(:,f1)<(nv(f1)/2) & tru_atomps(:,f2)>=(nv(f2)/2)  );
NR_21_id     = find(tru_atomps(:,f1)>=(nv(f1)/2) & tru_atomps(:,f2)<(nv(f2)/2) );
NR_22_id     = find(tru_atomps(:,f1)>=(nv(f1)/2) & tru_atomps(:,f2)>=(nv(f2)/2) );

atomps       = [atomps(NR_11_id,:);atomps(NR_12_id,:);atomps(NR_21_id,:);atomps(NR_22_id,:)];
```

0️⃣ Check the order of the orbitals ('orbitps') is sorted in the desired way throught the lines at the middle of the 'Supercell.m'. If there is no special needs for order of the atoms, comment these lines.

```Matlab
%%% ==================================================== %%%
%%% --- sort the orbitps (for Green ftn calculation) --- %%%
% orbitps      = sortrows(orbitps, 6);
frac_cord                   = 11:13;
frac_cord(ribbon_direction) = [];
%[~,index]                   = sortrows(floor(orbitps),frac_cord,'descend');
%orbitps                     = orbitps(index,:);

f1           = frac_cord(1);
f2           = frac_cord(2);
tru_orbitps  = orbitps;
tru_orbitps   = floor(tru_orbitps);

NR_11_id     = find(tru_orbitps(:,f1)<(nv(f1-6)/2) & tru_orbitps(:,f2)<(nv(f2-6)/2) );
NR_12_id     = find(tru_orbitps(:,f1)<(nv(f1-6)/2) & tru_orbitps(:,f2)>=(nv(f2-6)/2)  );
NR_21_id     = find(tru_orbitps(:,f1)>=(nv(f1-6)/2) & tru_orbitps(:,f2)<(nv(f2-6)/2) );
NR_22_id     = find(tru_orbitps(:,f1)>=(nv(f1-6)/2) & tru_orbitps(:,f2)>=(nv(f2-6)/2) );

orbitps      = [orbitps(NR_11_id,:);orbitps(NR_12_id,:);orbitps(NR_21_id,:);orbitps(NR_22_id,:)];
%%% ==================================================== %%%
```

0️⃣ If the super cell of a Sftn58sparse is desired, please check the normal of the slab is (001) in BR2D basis. The corresponding lines in 'SuperCell.m' are:

```Matlab
%% ------------------------------------------------------------
%%% if the supercell is built for slab with normal (001) in BR2D
if isSlab
    num_layer = max(ceil(ps(:,3)));
    ps(:,3)   = ps(:,3)/num_layer;
end
%% ------------------------------------------------------------
```

Such lines in 'SuperTBHmftn' are:

```Matlab
%% ------------------------------------------------------------
%%% if the supercell is built for slab with normal (001) in BR2D
if isSlab
    atomps(:,3)    = atomps(:,3)*num_layer;
    atomps(:,6:8)  = atomps(:,1:3)*BR;
    atomps(:,9:11) = atomps(:,6:8)/NBR;
    %orbitps(:,4:6) = atomps(:,9:11);
    orbitps(:,4:6) = ((orbitps(:,4:6)*NBR)/BR);
    orbitps(:,6)   = orbitps(:,6)*num_layer;
    orbitps(:,4:6) = (orbitps(:,4:6)*BR)/NBR;
end
%% ------------------------------------------------------------
```

1️⃣ Prepare the input file (input.txt) for the program "SuperCell.m"

```txt
SuperCell

wcal.vec1   = [2 0 0]                 % Transfromation matrix 
wcal.vec2   = [0 2 0]
wcal.vec3   = [0 0 1]
wcal.file   = ftn58sparse_cut_p.mat   % Input TB model 

endSuperCell
```

2️⃣ Run "SuperCell.m" to get the necessary infomation to construct the supercell TB model from the original TB model

```Matlab
>> SuperCell.m
```

After running this, a structure files "Unit.vasp" is generated. You can check "Unit.vasp" to see whether this is the supercell you want.

3️⃣ Run "SuperTBHmftn.m" to get the supercell TB model.

```Matlab
>> SuperTBHmftn.m
```

4️⃣ Check the band structure by running "band_ftn.m"

```Matlab
>> band_ftn.m
```

5️⃣ Change the "InFile" in "band_ftn.m" from "super_ftn58sparse.mat" to "ftn58sparse_cut_cdw.mat" for comparing the band structure w/o CDW and w/ CDW (supercell constructed from the prinstine monolayer TiSe2 ).

```Matlab
>> band_ftn.m
```

## 🔖 Tips

After running "SuperCell.m", a POSCAR file "Unit.vasp" will be generated. You can open it in VESTA to check that the geometry of the super cell, e.g. the Bravais vectors of the super cell, is the same as expected. Moreover, you can directly use it to do the first principle calculations for the super cell of this material.
