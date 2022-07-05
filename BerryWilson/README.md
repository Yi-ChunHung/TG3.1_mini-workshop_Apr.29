# 🧭 BerryWilson

This program calculate the Wilson loop for Berry connection of a real material countering along kx at fixed ky. The option for kz mesh is to observe the 3D TI (detailed principles please refer to the reference).

## 📓 [Reference](./Yu-2011-Equivalent-expression-of-z2-topological-invariant-using-non-Abelian-connection.pdf)

## 📓 [FurtherReading](./Fidkowski-2011-Model-characterization-of-gapless-edge-modes-of-TI-using-intermediate-BZ-functions.pdf)

## 📓 [See also](./2018%20Topological%20Crystalline%20Insulators.pdf)

## 🔰 Getting Started

0️⃣ Adjust the integrated direction in BerryWilson.m. If the ftn58sparse is 3D, also adjust the k-plane to calculate the Wilson loop in BerryWilson.m. The meaning of b1,b2,b3 are explained in (1.).

```Matlab
if isSlab 
    kvec = [b1(ii1),b2(ii2)]*2*pi;
else   
    kvec = [b1(ii1),b2(ii2),b3(ii3)]*2*pi;
end
```

1️⃣ Prepare the input file (input.txt) for the program "BerryWilson.m"

```txt
WilsonLoop

wcal.nVB    = 7                 % (nVB)th of occupied states 
wcal.ib1    = 40                % number of mesh point along the integrated direction on k-plane
wcal.ib2    = 40                % number of mesh point along another direction on k-plane
wcal.ib3    = 1                 % number of points along the rest direction to indicate the k-plane 
wcal.mode   = 0                 % mode=0: calculate loops on 2D EBZ ; mode=1: calculate loops on full 2D BZ
wcal.ref    = ftn58sparse.mat   % the BULK TB model
wcal.ncpu   = 6                 

endWilsonLoop
```

2️⃣ Run "BerryWilson.m" to calculate the phase of the Wilson loop for each ky

```Matlab
>> BerryWilson.m
```

3️⃣ Run "PlotLoops.m" to paint the Wilson loops for determination of Z2 and Chern number. (Please pay attention to the Effective Brilouin Zone for Z2 !!!)

```Matlab
>> PlotLoops.m
