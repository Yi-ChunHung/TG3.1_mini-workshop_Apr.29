# 🧭 STBHmftn

This program can create a finite slab tight-binding model with the (hkl) surface from the bulk tight-binding model.

## ❗ Warning

The format of Orbitps may greatly affect the validity of this program. Please pay attention on it for the ftn58sparse that may not comes from the standard process.

## 🔰 Getting Started

0️⃣ If you want to describe the plane through (intercept on Cartesian axis/normal of the plane in the basis of translation vectors), tunr (on/off) the 'isCar' option in 'SlabCell.m':

```Matlab
%% Check that [h k l] is in CARTESIAN representation (1) or in TRANSLATION representation (0)
%% (Added by YiChun 03/18/21)
isCAR = 1;
```

0️⃣ Check the order of the atoms ('superatomps') is sorted in the desired way through the lines at the bottom of 'SlabCell'. If there is no special needs for order of the atoms, comment these lines.

```Matlab
%%% -------------------------------------------------------- %%%
%%% --- sort the superatomps (for Green ftn calculation) --- %%%
%%% -------------------------------------------------------- %%%
frac_cord                   = 2;
superatomps                 = sortrows(superatomps,frac_cord,'ascend');
[~,index]                   = sortrows(floor(superatomps),frac_cord,'ascend');
superatomps(1:end)          = superatomps(index,:);
fprintf('sorted... \n')
```

1️⃣ Prepare the input file (input.txt) for the program "SlabCell.m"

```txt
SlabCell

wcal.nL     = 40                     % number of layers stacking along [hkl]
wcal.isCut  = 0                      % isCut = 1 to prepare the desired termination for the finite slab
wcal.LB     = [-0.1 0.5]             % if isCut = 1, LB is for the lower and upper boundaries. 
wcal.hkl    = [0 0 1]                % (hkl)
wcal.isConv = 0                      % For the case when the conventional unit vectors do not coincde with the Cartesian xyz
wcal.ref    = ftn58sparse_SnTe.mat   % name for the bulk TB model

endSlabCell
```

2️⃣ Run "SlabCell.m" to get the necessary infomation to construct the slab TB model from the bulk TB model

```Matlab
>> SlabCell.m
```

After running this, two structure files "Unit.vasp" and "Slab.vasp" are generated. You can check "Slab.vasp" to see whether this is the structure you want.

3️⃣ Run "STBHmftn.m" to get the slab TB model.

```Matlab
>> STBHmftn.m
```

4️⃣ Check the band structure by running "band_Sftn.m"

```Matlab
>> band_Sftn.m
```

## 📚 Notes

1️⃣ Script "MaxHop.m" will help to check the maximum hopping distance along the stakcing direction (the last index). To run the program "GreenSur.m", one will need the hopping information between layers that you can get from this program. Since the "STBHmftn.m" will reorient the surface normal as the third Bravis vector, one should focus on the third component of the results from "MaxHop.m".

2️⃣ If the atoms' postions are not as expected, please check the variable 'neighbor1', 'basis1', 'basis2' and 'basis3' in SlabCell.m. make sure that they are as expected. A +,- sign difference of them is not allowed in this program.

## 🔖 Tips

After running "SlabCell.m", a POSCAR file "Slab.vasp" will be generated. You can open it in VESTA to check that the geometry of the slab, e.g. the unit cell of the slab, is the same as expected. Moreover, you can directly use it to do the first principle calculations for the slab of this material.
