# 🧭 NRTBHmftn

Based on the program 'SuperCell', the program 'NRTBHmftn' generates a supercell TB model. Then, it truncate the hoppings along two specified directions to make the TB model for nano-ribbon.

## ❗ Caution

:one: Please pay attention to the unit cell of the nano-ribbon, especially when the direction that preserve translational symmetry does not proportion to any translation vector.

:two: The format of Orbitps will greatly affect the validity of this program. Please pay attention on it for the ftn58sparse that may not comes from the standard process.

## 🔰 Getting Started

1️⃣ Write the 'input.txt' for the information about the nano-ribbon you want to cut

```txt
SuperCell

wcal.vec1   = [1 1 0]
wcal.vec2   = [1 -1 0]
wcal.vec3   = [0 0 1]
wcal.file   = ftn58sparse.mat

wcal.ribdir = 3                       % the nth wcal.vec that preserves translational symmetry

endSuperCell
```

2️⃣ Run 'SuperCell.m' without any error occured to get the information about the chain supercell written in 'slab_info.mat'.

```Matlab
>> SuperCell.m
```

3️⃣ Run 'SuperTBHmftn.m' without any error occured to construct the TB model with the chain supercell written in 'super_ftn58sparse.mat'

```Matlab
>> SuperTBHmftn.m
```

4️⃣ Run 'NRTBHmftn.m' without any error occured. This program cuts off the hoppings on the edge of the TB model to get the TB model of chain written in 'Nano_Ribbon_ftn58sparse.mat'

```Matlab
>> NRTBHmftn.m
```

## 🔖 Tips

To find the translation vectors for the supercell to be cut into nano-ribbon, one may use the codes in ./FindTranslation :

1️⃣ Prepare the slab information of the side-surfaces on the nano-ribbon in "input_1" and "input_2", respectively.

```txt
SlabCell

wcal.nL     = 34
wcal.isCut  = 0
wcal.LB     = [-0.1 0.5]
wcal.hkl    = [1 1 4]
wcal.isConv = 0
wcal.ref    = ftn58sparse.mat

endSlabCell
```

2️⃣ Run the SlabCell_1.m and SlabCell_2.m to obtain the BR on the side-surfaces, respectively.

```Matlab
>> SlabCell.m
```

3️⃣ Run find_NRtranslation.m, the final results (i.e. the two translation vectors along the directions that would not preserve translational symmetry in nano-ribbon) will display in the command window.

```Matlab
>> find_NRtranslation.m
```

4️⃣ Eventually, after running "SuperCell.m", a POSCAR file "Unit.vasp" will be generated. You can open it in VESTA to check that the geometry of the nano-ribbon, e.g. the unit cell of the ribbon, is the same as expected. Moreover, you can directly use it to do the first principle calculations for the nano-ribbon of this material.
