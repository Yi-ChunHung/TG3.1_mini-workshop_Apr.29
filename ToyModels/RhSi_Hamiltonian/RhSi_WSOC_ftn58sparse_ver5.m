function ftn58sparse = RhSi_WSOC_ftn58sparse_ver5()
% ============================================
% -- translate RhSi_WSOC into ftn58sparse --
% ============================================

%% ========================================================================
ftn58sparse.System = 'RhSi_WSOC';
ftn58sparse.isToy  = 1; % no Ainfo, Orbitps, abc, BR, Nat
ftn58sparse.isSO   = 1;
ftn58sparse.norb   = 8;

%% ========================================================================

%--------------------------------------------------------------------------
Paulix = [0 1;1 0];
Pauliy = [0 -1i;1i 0];
Pauliz = [1 0;0 -1];
Ident  = [1 0;0 1];
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
C11 = kron(Ident, kron(Ident, Paulix));
C12 = kron(Ident, kron(Paulix, Paulix));
C13 = kron(Ident, kron(Paulix, Ident)); 
 
C21 = kron(Ident, kron(Pauliz, Pauliy));
C22 = kron(Ident, kron(Paulix, Pauliy));
C23 = kron(Ident, kron(Pauliy, Ident));
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
C3  = kron(kron(Ident,Ident),Ident);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
C_R11 = kron(Pauliy, kron(Pauliz, Pauliy));
C_R12 = kron(Pauliz, kron(Paulix, Pauliy));
C_R13 = kron(Paulix, kron(Pauliy, Ident));

C_R21 = kron(Pauliz, kron(Ident, Pauliy));
C_R22 = kron(Paulix, kron(Pauliy, Paulix));
C_R23 = kron(Pauliy, kron(Pauliy, Pauliz));

C_R31 = kron(Paulix, kron(Pauliz, Pauliy));
C_R32 = kron(Pauliy, kron(Paulix, Pauliy));
C_R33 = kron(Pauliz, kron(Pauliy, Ident));
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
C_S11 = kron(Paulix, kron(Ident, Paulix));
C_S12 = kron(Pauliy, kron(Paulix, Paulix));
C_S13 = kron(Pauliz, kron(Paulix, Ident));

C_S21 = kron(Pauliy, kron(Ident, Paulix));
C_S22 = kron(Pauliz, kron(Paulix, Paulix));
C_S23 = kron(Paulix, kron(Paulix, Ident));

C_S31 = kron(Pauliz, kron(Pauliz, Paulix));
C_S32 = kron(Paulix, kron(Pauliy, Pauliy));
C_S33 = kron(Pauliy, kron(Paulix, Pauliz));
%--------------------------------------------------------------------------

tvary = [0.55,-0.76,0.16,3.57354839937036,0,-0.03,0.01,-0.04,0,0];
norb  = 8;

%--------------------------------------------------------------------------    
% CxCy = cos(kp(1,1)/2)*cos(kp(1,2)/2);
% CyCz = cos(kp(1,2)/2)*cos(kp(1,3)/2);
% CzCx = cos(kp(1,1)/2)*cos(kp(1,3)/2);
    
% CxSy = cos(kp(1,1)/2)*sin(kp(1,2)/2);
% CySz = cos(kp(1,2)/2)*sin(kp(1,3)/2);
% CzSx = cos(kp(1,3)/2)*sin(kp(1,1)/2);
    
% SxSy = sin(kp(1,1)/2)*sin(kp(1,2)/2);
% SySz = sin(kp(1,2)/2)*sin(kp(1,3)/2);
% SzSx = sin(kp(1,1)/2)*sin(kp(1,3)/2);
    
% SxCy = sin(kp(1,1)/2)*cos(kp(1,2)/2);
% SyCz = sin(kp(1,2)/2)*cos(kp(1,3)/2);
% SzCx = sin(kp(1,3)/2)*cos(kp(1,1)/2);
%==========================================================================
% find ij indices for each term in Hamiltonian
% More specifically, the upper triangular part

[row_11,col_11] = find(C11);
row_11_temp     = row_11(row_11<col_11);
col_11     		= col_11(row_11<col_11);
row_11          = row_11_temp;

[row_12,col_12] = find(C12);
row_12_temp     = row_12(row_12<col_12);
col_12		    = col_12(row_12<col_12);
row_12          = row_12_temp;

[row_13,col_13] = find(C13);
row_13_temp     = row_13(row_13<col_13);
col_13          = col_13(row_13<col_13);
row_13          = row_13_temp;

[row_21,col_21] = find(C21);
row_21_temp          = row_21(row_21<col_21);
col_21        = col_21(row_21<col_21);
row_21          = row_21_temp;

[row_22,col_22] = find(C22);
row_22_temp          = row_22(row_22<col_22);
col_22          = col_22(row_22<col_22);
row_22          = row_22_temp;

[row_23,col_23] = find(C23);
row_23_temp          = row_23(row_23<col_23);
col_23          = col_23(row_23<col_23);
row_23          = row_23_temp;
 
 
[row_3,col_3]   = find(C3);   % on the diagonal of the 8x8
 
 
[row_R11,col_R11] = find(C_R11);
row_R11_temp          = row_R11(row_R11<col_R11);
col_R11          = col_R11(row_R11<col_R11);
row_R11          = row_R11_temp;

[row_R12,col_R12] = find(C_R12);
row_R12_temp          = row_R12(row_R12<col_R12);
col_R12          = col_R12(row_R12<col_R12);
row_R12          = row_R12_temp;

[row_R13,col_R13] = find(C_R13);
row_R13_temp          = row_R13(row_R13<col_R13);
col_R13          = col_R13(row_R13<col_R13);
row_R13          = row_R13_temp;

[row_R21,col_R21] = find(C_R21);
row_R21_temp          = row_R21(row_R21<col_R21);
col_R21          = col_R21(row_R21<col_R21);
row_R21          = row_R21_temp;

[row_R22,col_R22] = find(C_R22);
row_R22_temp          = row_R22(row_R22<col_R22);
col_R22          = col_R22(row_R22<col_R22);
row_R22          = row_R22_temp;

[row_R23,col_R23] = find(C_R23);
row_R23_temp          = row_R23(row_R23<col_R23);
col_R23         = col_R23(row_R23<col_R23);
row_R23          = row_R23_temp;

[row_R31,col_R31] = find(C_R31);
row_R31_temp          = row_R31(row_R31<col_R31);
col_R31          = col_R31(row_R31<col_R31);
row_R31          = row_R31_temp;

[row_R32,col_R32] = find(C_R32);
row_R32_temp          = row_R32(row_R32<col_R32);
col_R32          = col_R32(row_R32<col_R32);
row_R32          = row_R32_temp;

[row_R33,col_R33] = find(C_R33);
row_R33_temp          = row_R33(row_R33<col_R33);
col_R33          = col_R33(row_R33<col_R33);
row_R33          = row_R33_temp;

[row_S11,col_S11] = find(C_S11);
row_S11_temp          = row_S11(row_S11<col_S11);
col_S11          = col_S11(row_S11<col_S11);
row_S11          = row_S11_temp;

[row_S12,col_S12] = find(C_S12);
row_S12_temp          = row_S12(row_S12<col_S12);
col_S12          = col_S12(row_S12<col_S12);
row_S12          = row_S12_temp;

[row_S13,col_S13] = find(C_S13);
row_S13_temp          = row_S13(row_S13<col_S13);
col_S13          = col_S13(row_S13<col_S13);
row_S13          = row_S13_temp;

[row_S21,col_S21] = find(C_S21);
row_S21_temp          = row_S21(row_S21<col_S21);
col_S21          = col_S21(row_S21<col_S21);
row_S21          = row_S21_temp;

[row_S22,col_S22] = find(C_S22);
row_S22_temp          = row_S22(row_S22<col_S22);
col_S22          = col_S22(row_S22<col_S22);
row_S22          = row_S22_temp;

[row_S23,col_S23] = find(C_S23);
row_S23_temp          = row_S23(row_S23<col_S23);
col_S23          = col_S23(row_S23<col_S23);
row_S23          = row_S23_temp;

[row_S31,col_S31] = find(C_S31);
row_S31_temp          = row_S31(row_S31<col_S31);
col_S31          = col_S31(row_S31<col_S31);
row_S31          = row_S31_temp;

[row_S32,col_S32] = find(C_S32);
row_S32_temp          = row_S32(row_S32<col_S32);
col_S32          = col_S32(row_S32<col_S32);
row_S32          = row_S32_temp;

[row_S33,col_S33] = find(C_S33);
row_S33_temp          = row_S33(row_S33<col_S33);
col_S33          = col_S33(row_S33<col_S33);
row_S33          = row_S33_temp;
%--------------------------------------------------------------------------
% construct the hopping for each term in Hamiltonian

ftn58 = zeros(1,6);

% H1  = tvary(1)*(C11*CxCy+C12*CyCz+C13*CzCx);
dCxCy 	= [[0 0 0];[0 -1 0];[-1 0 0];[-1 -1 0]];
dCyCz 	= [[0 0 0];[0 -1 0];[0 0 -1];[0 -1 -1]];
dCxCz 	= [[0 0 0];[0 0 -1];[-1 0 0];[-1 0 -1]];

iijj_11  = [row_11,col_11];
iijjd_11 = combvec(iijj_11',dCxCy');
iijjd_11 = iijjd_11';
n        = size(iijjd_11,1);
ftn58	 = [ftn58;[iijjd_11(:,1:2) tvary(1)*(1/2)*(1/2)*ones(n,1) iijjd_11(:,3:5)]];

iijj_12  = [row_12,col_12];
iijjd_12 = combvec(iijj_12',dCyCz');
iijjd_12 = iijjd_12';
n        = size(iijjd_12,1);
ftn58	 = [ftn58;[iijjd_12(:,1:2) tvary(1)*(1/2)*(1/2)*ones(n,1) iijjd_12(:,3:5)]];

iijj_13  = [row_13,col_13];
iijjd_13 = combvec(iijj_13',dCxCz');
iijjd_13 = iijjd_13';
n        = size(iijjd_13,1);
ftn58	 = [ftn58;[iijjd_13(:,1:2) tvary(1)*(1/2)*(1/2)*ones(n,1) iijjd_13(:,3:5)]];

%H2  = tvary(2)*(C21*CxSy+C22*CySz+C23*CzSx);
dCxSy_p	= [[0 0 0];[-1 0 0]];
dCxSy_n = [[0 -1 0];[-1 -1 0]];
dCySz_p = [[0 0 0];[0 -1 0]];
dCySz_n = [[0 0 -1];[0 -1 -1]];
dCzSx_p = [[0 0 0];[0 0 -1]];
dCzSx_n = [[-1 0 0];[-1 0 -1]];

iijj_21     = [row_21,col_21];
iijjd_21_p  = combvec(iijj_21',dCxSy_p');
iijjd_21_p  = iijjd_21_p';
n        = size(iijjd_21_p,1);
ftn58	 = [ftn58;[iijjd_21_p(:,1:2) tvary(2)*(1/2)*(1/2i)*ones(n,1) iijjd_21_p(:,3:5)]];

iijj_21     = [row_21,col_21];
iijjd_21_n   = combvec(iijj_21',dCxSy_n');
iijjd_21_n  = iijjd_21_n';
n        = size(iijjd_21_n,1);
ftn58	 = [ftn58;[iijjd_21_n(:,1:2) (-1)*tvary(2)*(1/2)*(1/2i)*ones(n,1) iijjd_21_n(:,3:5)]];

iijj_22     = [row_22,col_22];
iijjd_22_p  = combvec(iijj_22',dCySz_p');
iijjd_22_p  = iijjd_22_p';
n        = size(iijjd_22_p,1);
ftn58	 = [ftn58;[iijjd_22_p(:,1:2) tvary(2)*(1/2)*(1/2i)*ones(n,1) iijjd_22_p(:,3:5)]];

iijj_22     = [row_22,col_22];
iijjd_22_n   = combvec(iijj_22',dCySz_n');
iijjd_22_n  = iijjd_22_n';
n        = size(iijjd_22_n,1);
ftn58	 = [ftn58;[iijjd_22_n(:,1:2) (-1)*tvary(2)*(1/2)*(1/2i)*ones(n,1) iijjd_22_n(:,3:5)]];

iijj_23     = [row_23,col_23];
iijjd_23_p  = combvec(iijj_23',dCzSx_p');
iijjd_23_p  = iijjd_23_p';
n        = size(iijjd_23_p,1);
ftn58	 = [ftn58;[iijjd_23_p(:,1:2) tvary(2)*(1/2)*(1/2i)*ones(n,1) iijjd_23_p(:,3:5)]];

iijj_23     = [row_23,col_23];
iijjd_23_n   = combvec(iijj_23',dCzSx_n');
iijjd_23_n  = iijjd_23_n';
n        = size(iijjd_23_n,1);
ftn58	 = [ftn58;[iijjd_23_n(:,1:2) (-1)*tvary(2)*(1/2)*(1/2i)*ones(n,1) iijjd_23_n(:,3:5)]];

%H3  = tvary(3)*C3*(cos(kp(1,1))+cos(kp(1,2))+cos(kp(1,3)));
dCx     = [[1 0 0];[-1 0 0]];
dCy     = [[0 1 0];[0 -1 0]];
dCz     = [[0 0 1];[0 0 -1]];
dCx_Cy_Cz = [dCx;dCy;dCz];

iijj_3  = [row_3,col_3];
iijjd_3 = combvec(iijj_3',dCx_Cy_Cz');
iijjd_3 = iijjd_3';
n        = size(iijjd_3,1);
ftn58	 = [ftn58;[iijjd_3(:,1:2) tvary(3)*(1/2)*ones(n,1) iijjd_3(:,3:5)]];

%VR1 = tvary(5)*(C_R11*CxCy+C_R12*CyCz+C_R13*CzCx);
iijj_R11  = [row_R11,col_R11];
iijjd_R11 = combvec(iijj_R11',dCxCy');
iijjd_R11 = iijjd_R11';
n        = size(iijjd_R11,1);
ftn58	 = [ftn58;[iijjd_R11(:,1:2) tvary(5)*(1/2)*(1/2)*ones(n,1) iijjd_R11(:,3:5)]];

iijj_R12  = [row_R12,col_R12];
iijjd_R12 = combvec(iijj_R12',dCyCz');
iijjd_R12 = iijjd_R12';
n        = size(iijjd_R12,1);
ftn58	 = [ftn58;[iijjd_R12(:,1:2) tvary(5)*(1/2)*(1/2)*ones(n,1) iijjd_R12(:,3:5)]];

iijj_R13  = [row_R13,col_R13];
iijjd_R13 = combvec(iijj_R13',dCxCz');
iijjd_R13 = iijjd_R13';
n        = size(iijjd_R13,1);
ftn58	 = [ftn58;[iijjd_R13(:,1:2) tvary(5)*(1/2)*(1/2)*ones(n,1) iijjd_R13(:,3:5)]];

%VR2 = tvary(6)*(C_R21*CxCy+C_R22*CyCz+C_R23*CzCx);
iijj_R21  = [row_R21,col_R21];
iijjd_R21 = combvec(iijj_R21',dCxCy');
iijjd_R21 = iijjd_R21';
n        = size(iijjd_R21,1);
ftn58	 = [ftn58;[iijjd_R21(:,1:2) tvary(6)*(1/2)*(1/2)*ones(n,1) iijjd_R21(:,3:5)]];

iijj_R22  = [row_R22,col_R22];
iijjd_R22 = combvec(iijj_R22',dCyCz');
iijjd_R22 = iijjd_R22';
n        = size(iijjd_R22,1);
ftn58	 = [ftn58;[iijjd_R22(:,1:2) tvary(6)*(1/2)*(1/2)*ones(n,1) iijjd_R22(:,3:5)]];

iijj_R23  = [row_R23,col_R23];
iijjd_R23 = combvec(iijj_R23',dCxCz');
iijjd_R23 = iijjd_R23';
n        = size(iijjd_R23,1);
ftn58	 = [ftn58;[iijjd_R23(:,1:2) tvary(6)*(1/2)*(1/2)*ones(n,1) iijjd_R23(:,3:5)]];

%VR3 = tvary(7)*(C_R31*SxSy+C_R32*SySz+C_R33*SzSx);
dSxSy_p	= [[0 0 0];[-1 -1 0]];
dSxSy_n = [[0 -1 0];[-1 0 0]];
dSySz_p	= [[0 0 0];[0 -1 -1]];
dSySz_n = [[0 -1 0];[0 0 -1]];
dSxSz_p	= [[0 0 0];[-1 0 -1]];
dSxSz_n = [[0 0 -1];[-1 0 0]];

iijj_R31    = [row_R31,col_R31];
iijjd_R31_p = combvec(iijj_R31',dSxSy_p');
iijjd_R31_p = iijjd_R31_p';
n        = size(iijjd_R31_p,1);
ftn58	 = [ftn58;[iijjd_R31_p(:,1:2) tvary(7)*(1/2i)*(1/2i)*ones(n,1) iijjd_R31_p(:,3:5)]];

iijj_R31    = [row_R31,col_R31];
iijjd_R31_n = combvec(iijj_R31',dSxSy_n');
iijjd_R31_n = iijjd_R31_n';
n        = size(iijjd_R31_n,1);
ftn58	 = [ftn58;[iijjd_R31_n(:,1:2) (-1)*tvary(7)*(1/2i)*(1/2i)*ones(n,1) iijjd_R31_n(:,3:5)]];

iijj_R32    = [row_R32,col_R32];
iijjd_R32_p = combvec(iijj_R32',dSySz_p');
iijjd_R32_p = iijjd_R32_p';
n        = size(iijjd_R32_p,1);
ftn58	 = [ftn58;[iijjd_R32_p(:,1:2) tvary(7)*(1/2i)*(1/2i)*ones(n,1) iijjd_R32_p(:,3:5)]];

iijj_R32    = [row_R32,col_R32];
iijjd_R32_n = combvec(iijj_R32',dSySz_n');
iijjd_R32_n = iijjd_R32_n';
n        = size(iijjd_R32_n,1);
ftn58	 = [ftn58;[iijjd_R32_n(:,1:2) (-1)*tvary(7)*(1/2i)*(1/2i)*ones(n,1) iijjd_R32_n(:,3:5)]];

iijj_R33    = [row_R33,col_R33];
iijjd_R33_p = combvec(iijj_R33',dSxSz_p');
iijjd_R33_p = iijjd_R33_p';
n        = size(iijjd_R33_p,1);
ftn58	 = [ftn58;[iijjd_R33_p(:,1:2) tvary(7)*(1/2i)*(1/2i)*ones(n,1) iijjd_R33_p(:,3:5)]];

iijj_R33    = [row_R33,col_R33];
iijjd_R33_n = combvec(iijj_R33',dSxSz_n');
iijjd_R33_n = iijjd_R33_n';
n        = size(iijjd_R33_n,1);
ftn58	 = [ftn58;[iijjd_R33_n(:,1:2) (-1)*tvary(7)*(1/2i)*(1/2i)*ones(n,1) iijjd_R33_n(:,3:5)]];

%VS1 = tvary(8)*(C_S11*SxCy+C_S12*SyCz+C_S13*SzCx);
dSxCy_p = [[0 0 0];[0 -1 0]];
dSxCy_n = [[-1 0 0];[-1 -1 0]];
dSyCz_p = [[0 0 0];[0 0 -1]];
dSyCz_n = [[0 -1 0];[0 -1 -1]];
dSzCx_p = [[0 0 0];[-1 0 0]];
dSzCx_n = [[0 0 -1];[-1 0 -1]];

iijj_S11    = [row_S11,col_S11];
iijjd_S11_p = combvec(iijj_S11',dSxCy_p');
iijjd_S11_p = iijjd_S11_p';
n        = size(iijjd_S11_p,1);
ftn58	 = [ftn58;[iijjd_S11_p(:,1:2) tvary(8)*(1/2i)*(1/2)*ones(n,1) iijjd_S11_p(:,3:5)]];

iijj_S11    = [row_S11,col_S11];
iijjd_S11_n = combvec(iijj_S11',dSxCy_n');
iijjd_S11_n = iijjd_S11_n';
n        = size(iijjd_S11_n,1);
ftn58	 = [ftn58;[iijjd_S11_n(:,1:2) (-1)*tvary(8)*(1/2i)*(1/2)*ones(n,1) iijjd_S11_n(:,3:5)]];

iijj_S12    = [row_S12,col_S12];
iijjd_S12_p = combvec(iijj_S12',dSyCz_p');
iijjd_S12_p = iijjd_S12_p';
n        = size(iijjd_S12_p,1);
ftn58	 = [ftn58;[iijjd_S12_p(:,1:2) tvary(8)*(1/2i)*(1/2)*ones(n,1) iijjd_S12_p(:,3:5)]];

iijj_S12    = [row_S12,col_S12];
iijjd_S12_n = combvec(iijj_S12',dSyCz_n');
iijjd_S12_n = iijjd_S12_n';
n        = size(iijjd_S12_n,1);
ftn58	 = [ftn58;[iijjd_S12_n(:,1:2) (-1)*tvary(8)*(1/2i)*(1/2)*ones(n,1) iijjd_S12_n(:,3:5)]];

iijj_S13    = [row_S13,col_S13];
iijjd_S13_p = combvec(iijj_S13',dSzCx_p');
iijjd_S13_p = iijjd_S13_p';
n        = size(iijjd_S13_p,1);
ftn58	 = [ftn58;[iijjd_S13_p(:,1:2) tvary(8)*(1/2i)*(1/2)*ones(n,1) iijjd_S13_p(:,3:5)]];

iijj_S13    = [row_S13,col_S13];
iijjd_S13_n = combvec(iijj_S13',dSzCx_n');
iijjd_S13_n = iijjd_S13_n';
n        = size(iijjd_S13_n,1);
ftn58	 = [ftn58;[iijjd_S13_n(:,1:2) (-1)*tvary(8)*(1/2i)*(1/2)*ones(n,1) iijjd_S13_n(:,3:5)]];

%VS2 = tvary(9)*(C_S21*CxSy+C_S22*CySz+C_S23*CzSx);
iijj_S21     = [row_S21,col_S21];
iijjd_S21_p  = combvec(iijj_S21',dCxSy_p');
iijjd_S21_p  = iijjd_S21_p';
n        = size(iijjd_S21_p,1);
ftn58	 = [ftn58;[iijjd_S21_p(:,1:2) tvary(9)*(1/2)*(1/2i)*ones(n,1) iijjd_S21_p(:,3:5)]];

iijj_S21     = [row_S21,col_S21];
iijjd_S21_n  = combvec(iijj_S21',dCxSy_n');
iijjd_S21_n  = iijjd_S21_n';
n        = size(iijjd_S21_n,1);
ftn58	 = [ftn58;[iijjd_S21_n(:,1:2) (-1)*tvary(9)*(1/2)*(1/2i)*ones(n,1) iijjd_S21_n(:,3:5)]];

iijj_S22     = [row_S22,col_S22];
iijjd_S22_p  = combvec(iijj_S22',dCySz_p');
iijjd_S22_p  = iijjd_S22_p';
n        = size(iijjd_S22_p,1);
ftn58	 = [ftn58;[iijjd_S22_p(:,1:2) tvary(9)*(1/2)*(1/2i)*ones(n,1) iijjd_S22_p(:,3:5)]];

iijj_S22     = [row_S22,col_S22];
iijjd_S22_n  = combvec(iijj_S22',dCySz_n');
iijjd_S22_n  = iijjd_S22_n';
n        = size(iijjd_S22_n,1);
ftn58	 = [ftn58;[iijjd_S22_n(:,1:2) (-1)*tvary(9)*(1/2)*(1/2i)*ones(n,1) iijjd_S22_n(:,3:5)]];

iijj_S23     = [row_S23,col_S23];
iijjd_S23_p  = combvec(iijj_S23',dCzSx_p');
iijjd_S23_p  = iijjd_S23_p';
n        = size(iijjd_S23_p,1);
ftn58	 = [ftn58;[iijjd_S23_p(:,1:2) tvary(9)*(1/2)*(1/2i)*ones(n,1) iijjd_S23_p(:,3:5)]];

iijj_S23     = [row_S23,col_S23];
iijjd_S23_n  = combvec(iijj_S23',dCzSx_n');
iijjd_S23_n  = iijjd_S23_n';
n        = size(iijjd_S23_n,1);
ftn58	 = [ftn58;[iijjd_S23_n(:,1:2) (-1)*tvary(9)*(1/2)*(1/2i)*ones(n,1) iijjd_S23_n(:,3:5)]];

%VS3 = tvary(10)*(C_S31*CxSy+C_S32*CySz+C_S33*CzSx);
iijj_S31     = [row_S31,col_S31];
iijjd_S31_p  = combvec(iijj_S31',dCxSy_p');
iijjd_S31_p  = iijjd_S31_p';
n        = size(iijjd_S31_p,1);
ftn58	 = [ftn58;[iijjd_S31_p(:,1:2) tvary(10)*(1/2)*(1/2i)*ones(n,1) iijjd_S31_p(:,3:5)]];

iijj_S31     = [row_S31,col_S31];
iijjd_S31_n  = combvec(iijj_S31',dCxSy_n');
iijjd_S31_n  = iijjd_S31_n';
n        = size(iijjd_S31_n,1);
ftn58	 = [ftn58;[iijjd_S31_n(:,1:2) (-1)*tvary(10)*(1/2)*(1/2i)*ones(n,1) iijjd_S31_n(:,3:5)]];

iijj_S32     = [row_S32,col_S32];
iijjd_S32_p  = combvec(iijj_S32',dCySz_p');
iijjd_S32_p  = iijjd_S32_p';
n        = size(iijjd_S32_p,1);
ftn58	 = [ftn58;[iijjd_S32_p(:,1:2) tvary(10)*(1/2)*(1/2i)*ones(n,1) iijjd_S32_p(:,3:5)]];

iijj_S32     = [row_S32,col_S32];
iijjd_S32_n  = combvec(iijj_S32',dCySz_n');
iijjd_S32_n  = iijjd_S32_n';
n        = size(iijjd_S32_n,1);
ftn58	 = [ftn58;[iijjd_S32_n(:,1:2) (-1)*tvary(10)*(1/2)*(1/2i)*ones(n,1) iijjd_S32_n(:,3:5)]];

iijj_S33     = [row_S33,col_S33];
iijjd_S33_p  = combvec(iijj_S33',dCzSx_p');
iijjd_S33_p  = iijjd_S33_p';
n        = size(iijjd_S33_p,1);
ftn58	 = [ftn58;[iijjd_S33_p(:,1:2) tvary(10)*(1/2)*(1/2i)*ones(n,1) iijjd_S33_p(:,3:5)]];

iijj_S33     = [row_S33,col_S33];
iijjd_S33_n  = combvec(iijj_S33',dCzSx_n');
iijjd_S33_n  = iijjd_S33_n';
n        = size(iijjd_S33_n,1);
ftn58	 = [ftn58;[iijjd_S33_n(:,1:2) (-1)*tvary(10)*(1/2)*(1/2i)*ones(n,1) iijjd_S33_n(:,3:5)]];

%==========================================================================
ij = ftn58(2:end,1:2);
tt = ftn58(2:end,3);
dd = ftn58(2:end,4:6);

% compensate the lower triangular part
v          = find(ij(:,1)==ij(:,2));
ij_temp    = ij;
tt_temp    = tt;
dd_temp    = dd;
ij_temp(v,:) = [];
tt_temp(v) = []; 
dd_temp(v,:) = [];
ftn58      = [ftn58;[flip(ij_temp,2) conj(tt_temp) (-1)*dd_temp]];

ftn58sparse.ij = ftn58(2:end,1:2);
ftn58sparse.tt = ftn58(2:end,3);
ftn58sparse.dd = ftn58(2:end,4:6);

save RhSi_WSOC_ftn58sparse_ver5.mat ftn58sparse

end