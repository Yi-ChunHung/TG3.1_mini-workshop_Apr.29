# 🧭 GreenSur_bulkRSRG

This program can calculate the surface states on a semi-inifite thin film by using the Green's function method with bulk information.

## 🔰 Getting Started

1️⃣ Prepare the input file in 'input.txt'.

```txt
GreenSur

wcal.Ni      = 20                     % number of iteration for RSRG 
wcal.ef      = 0
wcal.Emu     = [-1,1]                 % energy window to scan
wcal.Ne      = 200                    % # of sample points in energy window
wcal.np      = 6                      % # of cpus for parallel computation
wcal.nn      = 0.0001                 % broadening
wcal.kc      = 600                    % # of sample points in a k-path
wcal.kp1     = [0.0  0.0]             % vertices of the k-path
wcal.kp2     = [0.0  1.0]
wcal.surface = 1
wcal.hop_d   = 1
wcal.file    = ftn58sparse_vr0_g0.mat

endGreenSur
```

, in which surface (x=1, y=2, z=3) determines the surface plane and hop_d determines the termination of the surface.

```Matlab
%------------- if hop_d=2(semi-infinite at opposite surface)
%   BulkBulkBulk
%   BulkBulkBulk
%   BulkBulkBulk
%   BulkBulkBulk
%------------- if hop_d=1(semi-infinite at this surface)
```

2️⃣ Run "TBGreenSur_bulkRSRG.m" to get the Ek dispersion of the semi-infinite thin film.

```Matlab
>> TBGreenSur_bulkRSRG.m
```

3️⃣ Run "PlotEK.m" to plot the spectrum weight

```Matlab
>> PlotEk.m
```

## 🏁 [Example](./example)

The example contains the figure of overlaping the results calculated from two TB Hamiltonians, in which the code from the link below is used to do the trick.

> John Iversen (2022). freezeColors / unfreezeColors (https://www.mathworks.com/matlabcentral/fileexchange/7943-freezecolors-unfreezecolors), MATLAB Central File Exchange. Retrieved March 28, 2022.
