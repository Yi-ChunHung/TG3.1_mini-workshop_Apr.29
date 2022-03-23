# 🧭 BandPlot

The function bandplot_concise_ftn58sparse.m reads the file 'kpt.labelinfo.dat' for the k-path to plot the band structure of a given ftn58sparse. The default function is for a 3D system. The functions for lower dimensional system is indicated in their name.

## 🔰 Getting started

1️⃣ Edit the '@.labelinfo.dat' file. The '@' can be replaced by any words or characters. Same format is used also for lower dimensional system. Then, the function for 2D sytsem will read the first two column of k-points and the function for 1D sytem will read the first column of k-points. Same logic for bulk projection plots.

```txt
% label   sample point    k-points
Gamma     0               0    0   0    % The sample point must begin with 0
L         500             0.5  0   0
X         1000            0.5  0.5 0
Gamma     2000            0    0   0
```

2️⃣ run bandplot_concise_ftn58sparse.m

```Matlab
>> bandplot_concise_ftn58sparse(ftn58sparse, kpt_filename, PCD_orb_1,PCD_orb_2,PCD_orb_3,PCD_orb_4)
```

3️⃣ Load the desired 'Ekplot.mat' file in 'PlotDispersion.m'

```Matlab
load Ekplot.mat
```

4️⃣ Run 'PlotDispersion.m'

```Matlab
>> PlotDispersion(Ef,Name)
```

## 🚩 [band spintexture](./band%20spintexture)

The bandplot_concise_ftn58sparse_sp.m calculate the spintexture of a given component of spin on the band. Use 'PlotDispersionSP.m' to plot the results.

> - [ ] examples
> - [x] passed debug
> - [ ] exist future works 

## 🚩 [bulk projection](./bulk%20projection)

The bandplot_concise_ftn58sparse_ex_bp.m project the bulk band structure onto lower dimensional BZ. The code assume the lower dimensional BZ is to describe the **broken periodicities along reciprocal lattice vectors**.

> - [x] examples
> - [x] passed debug
> - [ ] exist future works

## 🚩 [interface band](./interface%20band)

The bandplot_concise_ftn58sparse_interface.m provide a module for constructing a special kind of interface and plotting its band structure.

> - [ ] examples
> - [x] passed debug
> - [ ] exist future works

## 🏁 [Example](./example)

This is an example of 3D material by file bandplot_concise_ftn58sparse.m. In example, the material GaAs are constructed.
