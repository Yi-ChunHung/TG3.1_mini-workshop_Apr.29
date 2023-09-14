function ftn58sparse = RhSi_WSOC_ftn58sparse()
% ============================================
% -- translate RhSi_WSOC into ftn58sparse --
% ============================================

%% ========================================================================
ftn58sparse.System = 'RhSi_WSOC';
ftn58sparse.isToy  = 1; % no Ainfo, Orbitps, abc, BR, Nat
ftn58sparse.isSO   = 1;
ftn58sparse.norb   = 8;

[i,j]   = meshgrid(1:8,1:8);
ii      = reshape(i,[],1);
jj      = reshape(j,[],1);
iijj 	= [ii,jj];
 
d       = [[0 0 0];[0 -1 0];[-1 0 0];[-1 -1 0];...
		   [0 -1 0];[0 0 -1];[0 -1 -1];...
		   [0 0 -1];[-1 0 0];[-1 0 -1]];

iijjd 	= combvec(iijj',d');
iijjd   = iijjd';
ij      = iijjd(:,1:2);
dd      = iijjd(:,3:5); 
nbonds  = length(ij(:,1));

ftn58sparse.ij = ij;
ftn58sparse.dd = dd;
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
norb  = ftn58sparse.norb;

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

[row_11,col_11] = find(C11);
[row_12,col_12] = find(C12);
[row_13,col_13] = find(C13);
[row_21,col_21] = find(C21);
[row_22,col_22] = find(C22);
[row_23,col_23] = find(C23);
 
[row_3,col_3]   = find(C3);
 
[row_R11,col_R11] = find(C_R11);
[row_R12,col_R12] = find(C_R12);
[row_R13,col_R13] = find(C_R13);
[row_R21,col_R21] = find(C_R21);
[row_R22,col_R22] = find(C_R22);
[row_R23,col_R23] = find(C_R23);
[row_R31,col_R31] = find(C_R31);
[row_R32,col_R32] = find(C_R32);
[row_R33,col_R33] = find(C_R33);
 
[row_S11,col_S11] = find(C_S11);
[row_S12,col_S12] = find(C_S12);
[row_S13,col_S13] = find(C_S13);
[row_S21,col_S21] = find(C_S21);
[row_S22,col_S22] = find(C_S22);
[row_S23,col_S23] = find(C_S23);
[row_S31,col_S31] = find(C_S31);
[row_S32,col_S32] = find(C_S32);
[row_S33,col_S33] = find(C_S33);

%--------------------------------------------------------------------------
% construct the hopping for each term in Hamiltonian

tt = zeros(nbonds,1);

% H1  = tvary(1)*(C11*CxCy+C12*CyCz+C13*CzCx);
dCxCy 	= [[0 0 0];[0 -1 0];[-1 0 0];[-1 -1 0]];
dCyCz 	= [[0 0 0];[0 -1 0];[0 0 -1];[0 -1 -1]];
dCxCz 	= [[0 0 0];[0 0 -1];[-1 0 0];[-1 0 -1]];

iijj_11  = [row_11,col_11];
iijjd_11 = combvec(iijj_11',dCxCy');
iijjd_11 = iijjd_11';
v_11	 = ismember(iijjd,iijjd_11,'rows'); 
w_11	 = find(v_11);
tt(w_11) = tt(w_11) + tvary(1)*(1/2)*(1/2);

iijj_12  = [row_12,col_12];
iijjd_12 = combvec(iijj_12',dCyCz');
iijjd_12 = iijjd_12';
v_12	 = ismember(iijjd,iijjd_12,'rows'); 
w_12	 = find(v_12);
tt(w_12) = tt(w_12) + tvary(1)*(1/2)*(1/2);

iijj_13  = [row_13,col_13];
iijjd_13 = combvec(iijj_13',dCxCz');
iijjd_13 = iijjd_13';
v_13	 = ismember(iijjd,iijjd_13,'rows'); 
w_13	 = find(v_13);
tt(w_13) = tt(w_13) + tvary(1)*(1/2)*(1/2);

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
v_21_p		= ismember(iijjd,iijjd_21_p,'rows'); 
w_21_p	 	= find(v_21_p);
tt(w_21_p) 	= tt(w_21_p) + tvary(2)*(1/2)*(1/(2i));

iijj_21     = [row_21,col_21];
iijjd_21_n   = combvec(iijj_21',dCxSy_n');
iijjd_21_n  = iijjd_21_n';
v_21_n	 	= ismember(iijjd,iijjd_21_n,'rows'); 
w_21_n	 	= find(v_21_n);
tt(w_21_n) 	= tt(w_21_n) + (-1)*tvary(2)*(1/2)*(1/(2i));

iijj_22     = [row_22,col_22];
iijjd_22_p  = combvec(iijj_22',dCySz_p');
iijjd_22_p  = iijjd_22_p';
v_22_p	 	= ismember(iijjd,iijjd_22_p,'rows'); 
w_22_p	 	= find(v_22_p);
tt(w_22_p) 	= tt(w_22_p) + tvary(2)*(1/2)*(1/(2i));

iijj_22     = [row_22,col_22];
iijjd_22_n   = combvec(iijj_22',dCySz_n');
iijjd_22_n  = iijjd_22_n';
v_22_n	 	= ismember(iijjd,iijjd_22_n,'rows'); 
w_22_n	 	= find(v_22_n);
tt(w_22_n) 	= tt(w_22_n) + (-1)*tvary(2)*(1/2)*(1/(2i));

iijj_23     = [row_23,col_23];
iijjd_23_p  = combvec(iijj_23',dCzSx_p');
iijjd_23_p  = iijjd_23_p';
v_23_p	 	= ismember(iijjd,iijjd_23_p,'rows'); 
w_23_p	 	= find(v_23_p);
tt(w_23_p) 	= tt(w_23_p) + tvary(2)*(1/2)*(1/(2i));

iijj_23     = [row_23,col_23];
iijjd_23_n   = combvec(iijj_23',dCzSx_n');
iijjd_23_n  = iijjd_23_n';
v_23_n	 	= ismember(iijjd,iijjd_23_n,'rows'); 
w_23_n	 	= find(v_23_n);
tt(w_23_n) 	= tt(w_23_n) + (-1)*tvary(2)*(1/2)*(1/(2i));

%H3  = tvary(3)*C3*(cos(kp(1,1))+cos(kp(1,2))+cos(kp(1,3)));
dCx     = [[1 0 0];[-1 0 0]];
dCy     = [[0 1 0];[0 -1 0]];
dCz     = [[0 0 1];[0 0 -1]];
dCx_Cy_Cz = [dCx;dCy;dCz];

iijj_3  = [row_3,col_3];
iijjd_3 = combvec(iijj_3',dCx_Cy_Cz');
iijjd_3 = iijjd_3';
v_3	    = ismember(iijjd,iijjd_3,'rows'); 
w_3	    = find(v_3);
tt(w_3) = tt(w_3) + tvary(3)*(1/2);

%VR1 = tvary(5)*(C_R11*CxCy+C_R12*CyCz+C_R13*CzCx);
iijj_R11  = [row_R11,col_R11];
iijjd_R11 = combvec(iijj_R11',dCxCy');
iijjd_R11 = iijjd_R11';
v_R11	 = ismember(iijjd,iijjd_R11,'rows'); 
w_R11	 = find(v_R11);
tt(w_R11) = tt(w_R11) + tvary(5)*(1/2)*(1/2);

iijj_R12  = [row_R12,col_R12];
iijjd_R12 = combvec(iijj_R12',dCyCz');
iijjd_R12 = iijjd_R12';
v_R12	 = ismember(iijjd,iijjd_R12,'rows'); 
w_R12	 = find(v_R12);
tt(w_R12) = tt(w_R12) + tvary(5)*(1/2)*(1/2);

iijj_R13  = [row_R13,col_R13];
iijjd_R13 = combvec(iijj_R13',dCxCz');
iijjd_R13 = iijjd_R13';
v_R13	 = ismember(iijjd,iijjd_R13,'rows'); 
w_R13	 = find(v_R13);
tt(w_R13) = tt(w_R13) + tvary(5)*(1/2)*(1/2);

%VR2 = tvary(6)*(C_R21*CxCy+C_R22*CyCz+C_R23*CzCx);
iijj_R21  = [row_R21,col_R21];
iijjd_R21 = combvec(iijj_R21',dCxCy');
iijjd_R21 = iijjd_R21';
v_R21	 = ismember(iijjd,iijjd_R21,'rows'); 
w_R21	 = find(v_R21);
tt(w_R21) = tt(w_R21) + tvary(6)*(1/2)*(1/2);

iijj_R22  = [row_R22,col_R22];
iijjd_R22 = combvec(iijj_R22',dCyCz');
iijjd_R22 = iijjd_R22';
v_R22	 = ismember(iijjd,iijjd_R22,'rows'); 
w_R22	 = find(v_R22);
tt(w_R22) = tt(w_R22) + tvary(6)*(1/2)*(1/2);

iijj_R23  = [row_R23,col_R23];
iijjd_R23 = combvec(iijj_R23',dCxCz');
iijjd_R23 = iijjd_R23';
v_R23	 = ismember(iijjd,iijjd_R23,'rows'); 
w_R23	 = find(v_R23);
tt(w_R23) = tt(w_R23) + tvary(6)*(1/2)*(1/2);

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
v_R31_p	    = ismember(iijjd,iijjd_R31_p,'rows'); 
w_R31_p	    = find(v_R31_p);
tt(w_R31_p) = tt(w_R31_p) + tvary(7)*(1/2i)*(1/2i);

iijj_R31    = [row_R31,col_R31];
iijjd_R31_n = combvec(iijj_R31',dSxSy_n');
iijjd_R31_n = iijjd_R31_n';
v_R31_n	 	= ismember(iijjd,iijjd_R31_n,'rows'); 
w_R31_n	 	= find(v_R31_n);
tt(w_R31_n)	= tt(w_R31_n) + (-1)*tvary(7)*(1/2i)*(1/2i);

iijj_R32    = [row_R32,col_R32];
iijjd_R32_p = combvec(iijj_R32',dSySz_p');
iijjd_R32_p = iijjd_R32_p';
v_R32_p	    = ismember(iijjd,iijjd_R32_p,'rows'); 
w_R32_p	    = find(v_R32_p);
tt(w_R32_p) = tt(w_R32_p) + tvary(7)*(1/2i)*(1/2i);

iijj_R32    = [row_R32,col_R32];
iijjd_R32_n = combvec(iijj_R32',dSySz_n');
iijjd_R32_n = iijjd_R32_n';
v_R32_n	 	= ismember(iijjd,iijjd_R32_n,'rows'); 
w_R32_n	 	= find(v_R32_n);
tt(w_R32_n)	= tt(w_R32_n) + (-1)*tvary(7)*(1/2i)*(1/2i);

iijj_R33    = [row_R33,col_R33];
iijjd_R33_p = combvec(iijj_R33',dSxSz_p');
iijjd_R33_p = iijjd_R33_p';
v_R33_p	    = ismember(iijjd,iijjd_R33_p,'rows'); 
w_R33_p	    = find(v_R33_p);
tt(w_R33_p) = tt(w_R33_p) + tvary(7)*(1/2i)*(1/2i);

iijj_R33    = [row_R33,col_R33];
iijjd_R33_n = combvec(iijj_R33',dSxSz_n');
iijjd_R33_n = iijjd_R33_n';
v_R33_n	 	= ismember(iijjd,iijjd_R33_n,'rows'); 
w_R33_n	 	= find(v_R33_n);
tt(w_R33_n)	= tt(w_R33_n) + (-1)*tvary(7)*(1/2i)*(1/2i);

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
v_S11_p	    = ismember(iijjd,iijjd_S11_p,'rows'); 
w_S11_p	    = find(v_S11_p);
tt(w_S11_p) = tt(w_S11_p) + tvary(8)*(1/2i)*(1/2);

iijj_S11    = [row_S11,col_S11];
iijjd_S11_n = combvec(iijj_S11',dSxCy_n');
iijjd_S11_n = iijjd_S11_n';
v_S11_n	 	= ismember(iijjd,iijjd_S11_n,'rows'); 
w_S11_n	 	= find(v_S11_n);
tt(w_S11_n)	= tt(w_S11_n) + (-1)*tvary(8)*(1/2i)*(1/2);

iijj_S12    = [row_S12,col_S12];
iijjd_S12_p = combvec(iijj_S12',dSyCz_p');
iijjd_S12_p = iijjd_S12_p';
v_S12_p	    = ismember(iijjd,iijjd_S12_p,'rows'); 
w_S12_p	    = find(v_S12_p);
tt(w_S12_p) = tt(w_S12_p) + tvary(8)*(1/2i)*(1/2);

iijj_S12    = [row_S12,col_S12];
iijjd_S12_n = combvec(iijj_S12',dSyCz_n');
iijjd_S12_n = iijjd_S12_n';
v_S12_n	 	= ismember(iijjd,iijjd_S12_n,'rows'); 
w_S12_n	 	= find(v_S12_n);
tt(w_S12_n)	= tt(w_S12_n) + (-1)*tvary(8)*(1/2i)*(1/2);

iijj_S13    = [row_S13,col_S13];
iijjd_S13_p = combvec(iijj_S13',dSzCx_p');
iijjd_S13_p = iijjd_S13_p';
v_S13_p	    = ismember(iijjd,iijjd_S13_p,'rows'); 
w_S13_p	    = find(v_S13_p);
tt(w_S13_p) = tt(w_S13_p) + tvary(8)*(1/2i)*(1/2);

iijj_S13    = [row_S13,col_S13];
iijjd_S13_n = combvec(iijj_S13',dSzCx_n');
iijjd_S13_n = iijjd_S13_n';
v_S13_n	 	= ismember(iijjd,iijjd_S13_n,'rows'); 
w_S13_n	 	= find(v_S13_n);
tt(w_S13_n)	= tt(w_S13_n) + (-1)*tvary(8)*(1/2i)*(1/2);

%VS2 = tvary(9)*(C_S21*CxSy+C_S22*CySz+C_S23*CzSx);
iijj_S21     = [row_S21,col_S21];
iijjd_S21_p  = combvec(iijj_S21',dCxSy_p');
iijjd_S21_p  = iijjd_S21_p';
v_S21_p		 = ismember(iijjd,iijjd_S21_p,'rows'); 
w_S21_p	 	 = find(v_S21_p);
tt(w_S21_p)  = tt(w_S21_p) + tvary(9)*(1/2)*(1/(2i));

iijj_S21     = [row_S21,col_S21];
iijjd_S21_n  = combvec(iijj_S21',dCxSy_n');
iijjd_S21_n  = iijjd_S21_n';
v_S21_n	 	 = ismember(iijjd,iijjd_S21_n,'rows'); 
w_S21_n      = find(v_S21_n);
tt(w_S21_n)  = tt(w_S21_n) + (-1)*tvary(9)*(1/2)*(1/(2i));

iijj_S22     = [row_S22,col_S22];
iijjd_S22_p  = combvec(iijj_S22',dCySz_p');
iijjd_S22_p  = iijjd_S22_p';
v_S22_p	 	 = ismember(iijjd,iijjd_S22_p,'rows'); 
w_S22_p	 	 = find(v_S22_p);
tt(w_S22_p)  = tt(w_S22_p) + tvary(9)*(1/2)*(1/(2i));

iijj_S22     = [row_S22,col_S22];
iijjd_S22_n  = combvec(iijj_S22',dCySz_n');
iijjd_S22_n  = iijjd_S22_n';
v_S22_n	 	 = ismember(iijjd,iijjd_S22_n,'rows'); 
w_S22_n	 	 = find(v_S22_n);
tt(w_S22_n)  = tt(w_S22_n) + (-1)*tvary(9)*(1/2)*(1/(2i));

iijj_S23     = [row_S23,col_S23];
iijjd_S23_p  = combvec(iijj_S23',dCzSx_p');
iijjd_S23_p  = iijjd_S23_p';
v_S23_p      = ismember(iijjd,iijjd_S23_p,'rows'); 
w_S23_p	 	 = find(v_S23_p);
tt(w_S23_p)  = tt(w_S23_p) + tvary(9)*(1/2)*(1/(2i));

iijj_S23     = [row_S23,col_S23];
iijjd_S23_n  = combvec(iijj_S23',dCzSx_n');
iijjd_S23_n  = iijjd_S23_n';
v_S23_n	 	 = ismember(iijjd,iijjd_S23_n,'rows'); 
w_S23_n	 	 = find(v_S23_n);
tt(w_S23_n)  = tt(w_S23_n) + (-1)*tvary(9)*(1/2)*(1/(2i));

%VS3 = tvary(10)*(C_S31*CxSy+C_S32*CySz+C_S33*CzSx);
iijj_S31     = [row_S31,col_S31];
iijjd_S31_p  = combvec(iijj_S31',dCxSy_p');
iijjd_S31_p  = iijjd_S31_p';
v_S31_p		 = ismember(iijjd,iijjd_S31_p,'rows'); 
w_S31_p	 	 = find(v_S31_p);
tt(w_S31_p)  = tt(w_S31_p) + tvary(10)*(1/2)*(1/(2i));

iijj_S31     = [row_S31,col_S31];
iijjd_S31_n  = combvec(iijj_S31',dCxSy_n');
iijjd_S31_n  = iijjd_S31_n';
v_S31_n	 	 = ismember(iijjd,iijjd_S31_n,'rows'); 
w_S31_n	 	 = find(v_S31_n);
tt(w_S31_n)  = tt(w_S31_n) + (-1)*tvary(10)*(1/2)*(1/(2i));

iijj_S32     = [row_S32,col_S32];
iijjd_S32_p  = combvec(iijj_S32',dCySz_p');
iijjd_S32_p  = iijjd_S32_p';
v_S32_p	 	 = ismember(iijjd,iijjd_S32_p,'rows'); 
w_S32_p	 	 = find(v_S32_p);
tt(w_S32_p)  = tt(w_S32_p) + tvary(10)*(1/2)*(1/(2i));

iijj_S32     = [row_S32,col_S32];
iijjd_S32_n  = combvec(iijj_S32',dCySz_n');
iijjd_S32_n  = iijjd_S32_n';
v_S32_n	 	 = ismember(iijjd,iijjd_S32_n,'rows'); 
w_S32_n	 	 = find(v_S32_n);
tt(w_S32_n)  = tt(w_S32_n) + (-1)*tvary(10)*(1/2)*(1/(2i));

iijj_S33     = [row_S33,col_S33];
iijjd_S33_p  = combvec(iijj_S33',dCzSx_p');
iijjd_S33_p  = iijjd_S33_p';
v_S33_p	 	 = ismember(iijjd,iijjd_S33_p,'rows'); 
w_S33_p	 	 = find(v_S33_p);
tt(w_S33_p)  = tt(w_S33_p) + tvary(10)*(1/2)*(1/(2i));

iijj_S33     = [row_S33,col_S33];
iijjd_S33_n  = combvec(iijj_S33',dCzSx_n');
iijjd_S33_n  = iijjd_S33_n';
v_S33_n	 	 = ismember(iijjd,iijjd_S33_n,'rows'); 
w_S33_n	 	 = find(v_S33_n);
tt(w_S33_n)  = tt(w_S33_n) + (-1)*tvary(10)*(1/2)*(1/(2i));


%==========================================================================
ftn58sparse.tt = tt;

save RhSi_WSOC_ftn58sparse_ver2.mat ftn58sparse

end