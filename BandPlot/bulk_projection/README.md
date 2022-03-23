# 🧭 bandplot_concise_ftn58sparse_ex_bp.m

The function reads the file 'kpt.labelinfo.dat' for the k-path to plot the projected band structure of a given ftn58sparse.

- _**bandplot_concise_ftn58sparse_ex_bp.m**_:  
  project a 3D bulk band structure to a 2D BZ
  
- _**bandplot_concise_ftn58sparse2D_ex_bp.m**_:  
  project a 2D bulk band structure to a 1D BZ
  
- _**bandplot_concise_ftn58sparse3Dto1D_ex_bp.m**_:  
  project a 3D bulk band structure to a 1D BZ

## 🔰 Getting started

0️⃣ Check the directions of the projection  (or the direction of the periodic direction) is the desired one:

```Matlab
bandplot_concise_ftn58sparse_ex_bp(ftn58sparse, Ef, isSP, "BPdirection" , PCD_orb_1,PCD_orb_2,PCD_orb_3,PCD_orb_4,Name)
```

```Matlab
bandplot_concise_ftn58sparse3Dto1D_ex_bp(ftn58sparse, Ef, isSP, "Pdirection", PCD_orb_1,PCD_orb_2,PCD_orb_3,PCD_orb_4,Name)
```

Note that 'BPdirection' and 'Pdirection' = 1, 2, 3 are numbers.

1️⃣ Edit the '@.labelinfo.dat' file. The '@' can be replaced by any words or characters. Same format is used also for lower dimensional system. Then, the function for 2D sytsem will read the first two column of k-points and the function for 1D sytem will read the first column of k-points. Same logic for bulk projection plots.

```txt
% label   sample point    k-points
Gamma     0               0    0   0    % The sample point must begin with 0
X         500             0.5  0   0
L         1000            0.5  0.5 0
Gamma     2000            0    0   0
```

2️⃣ Load the desired 'Ekplot.mat' file in 'PlotDispersion.m'

```Matlab
load Ekplot.mat
```

3️⃣ Run 'PlotDispersion.m' ( or 'PlotDispersion(3Dto1D).m')

```Matlab
>> PlotDispersion(Ef,Name)
```

## 🏁 [Example](./example)

This is an example of 3D material project to 2D surface by file bandplot_concise_ftn58sparse_ex_bp.m. In example, the material GaAs are constructed.
