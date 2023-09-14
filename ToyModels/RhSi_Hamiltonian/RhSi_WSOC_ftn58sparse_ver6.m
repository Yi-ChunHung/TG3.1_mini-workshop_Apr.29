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
% construct the hopping for each term in Hamiltonian

ftn58 = zeros(1,6);

% H1  = tvary(1)*(C11*CxCy+C12*CyCz+C13*CzCx);
dCxCy 	= [[0 0 0];[0 -1 0];[-1 0 0];[-1 -1 0]];
dCyCz 	= [[0 0 0];[0 -1 0];[0 0 -1];[0 -1 -1]];
dCzCx 	= [[0 0 0];[0 0 -1];[-1 0 0];[-1 0 -1]];

[ii_11,jj_11,tt_11] = find(sparse(C11));
tt_11               = (1/2)*(1/2)*tvary(1)*tt_11;
iijjtt_11           = [ii_11,jj_11,tt_11];
iijjttdd_11         = combvec(iijjtt_11',dCxCy');
ftn58               = [ftn58;iijjttdd_11'];

[ii_12,jj_12,tt_12] = find(sparse(C12));
tt_12               = (1/2)*(1/2)*tvary(1)*tt_12;
iijjtt_12           = [ii_12,jj_12,tt_12];
iijjttdd_12         = combvec(iijjtt_12',dCyCz');
ftn58               = [ftn58;iijjttdd_12'];

[ii_13,jj_13,tt_13] = find(sparse(C13));
tt_13               = (1/2)*(1/2)*tvary(1)*tt_13;
iijjtt_13           = [ii_13,jj_13,tt_13];
iijjttdd_13         = combvec(iijjtt_13',dCzCx');
ftn58               = [ftn58;iijjttdd_13'];

%H2  = tvary(2)*(C21*CxSy+C22*CySz+C23*CzSx);
dCxSy_p	= [[0 0 0];[-1 0 0]];
dCxSy_n = [[0 -1 0];[-1 -1 0]];
dCySz_p = [[0 0 0];[0 -1 0]];
dCySz_n = [[0 0 -1];[0 -1 -1]];
dCzSx_p = [[0 0 0];[0 0 -1]];
dCzSx_n = [[-1 0 0];[-1 0 -1]];

[ii_21_p,jj_21_p,tt_21_p] = find(sparse(C21));
tt_21_p                   = (1/2)*(1/2i)*tvary(2)*tt_21_p;
iijjtt_21_p               = [ii_21_p,jj_21_p,tt_21_p];
iijjttdd_21_p     		  = combvec(iijjtt_21_p',dCxSy_p');
ftn58               	  = [ftn58;iijjttdd_21_p'];

[ii_21_n,jj_21_n,tt_21_n] = find(sparse(C21));
tt_21_n                   = (-1)*(1/2)*(1/2i)*tvary(2)*tt_21_n;
iijjtt_21_n               = [ii_21_n,jj_21_n,tt_21_n];
iijjttdd_21_n     		  = combvec(iijjtt_21_n',dCxSy_n');
ftn58               	  = [ftn58;iijjttdd_21_n'];

[ii_22_p,jj_22_p,tt_22_p] = find(sparse(C22));
tt_22_p                   = (1/2)*(1/2i)*tvary(2)*tt_22_p;
iijjtt_22_p               = [ii_22_p,jj_22_p,tt_22_p];
iijjttdd_22_p     		  = combvec(iijjtt_22_p',dCySz_p');
ftn58               	  = [ftn58;iijjttdd_22_p'];

[ii_22_n,jj_22_n,tt_22_n] = find(sparse(C22));
tt_22_n                   = (-1)*(1/2)*(1/2i)*tvary(2)*tt_22_n;
iijjtt_22_n               = [ii_22_n,jj_22_n,tt_22_n];
iijjttdd_22_n     		  = combvec(iijjtt_22_n',dCySz_n');
ftn58               	  = [ftn58;iijjttdd_22_n'];

[ii_23_p,jj_23_p,tt_23_p] = find(sparse(C23));
tt_23_p                   = (1/2)*(1/2i)*tvary(2)*tt_23_p;
iijjtt_23_p               = [ii_23_p,jj_23_p,tt_23_p];
iijjttdd_23_p     		  = combvec(iijjtt_23_p',dCzSx_p');
ftn58               	  = [ftn58;iijjttdd_23_p'];

[ii_23_n,jj_23_n,tt_23_n] = find(sparse(C23));
tt_23_n                   = (-1)*(1/2)*(1/2i)*tvary(2)*tt_23_n;
iijjtt_23_n               = [ii_23_n,jj_23_n,tt_23_n];
iijjttdd_23_n     		  = combvec(iijjtt_23_n',dCzSx_n');
ftn58               	  = [ftn58;iijjttdd_23_n'];

%H3  = tvary(3)*C3*(cos(kp(1,1))+cos(kp(1,2))+cos(kp(1,3)));
dCx     = [[1 0 0];[-1 0 0]];
dCy     = [[0 1 0];[0 -1 0]];
dCz     = [[0 0 1];[0 0 -1]];
dCx_Cy_Cz = [dCx;dCy;dCz];

[ii_3,jj_3,tt_3] = find(sparse(C3));
tt_3             = (1/2)*tvary(3)*tt_3;
iijjtt_3         = [ii_3,jj_3,tt_3];
iijjttdd_3       = combvec(iijjtt_3',dCx_Cy_Cz');
ftn58            = [ftn58;iijjttdd_3'];

%VR1 = tvary(5)*(C_R11*CxCy+C_R12*CyCz+C_R13*CzCx);
[ii_R11,jj_R11,tt_R11] = find(sparse(C_R11));
tt_R11                 = (1/2)*(1/2)*tvary(5)*tt_R11;
iijjtt_R11             = [ii_R11,jj_R11,tt_R11];
iijjttdd_R11           = combvec(iijjtt_R11',dCxCy');
ftn58                  = [ftn58;iijjttdd_R11'];

[ii_R12,jj_R12,tt_R12] = find(sparse(C_R12));
tt_R12                 = (1/2)*(1/2)*tvary(5)*tt_R12;
iijjtt_R12             = [ii_R12,jj_R12,tt_R12];
iijjttdd_R12           = combvec(iijjtt_R12',dCyCz');
ftn58                  = [ftn58;iijjttdd_R12'];

[ii_R13,jj_R13,tt_R13] = find(sparse(C_R13));
tt_R13                 = (1/2)*(1/2)*tvary(5)*tt_R13;
iijjtt_R13             = [ii_R13,jj_R13,tt_R13];
iijjttdd_R13           = combvec(iijjtt_R13',dCzCx');
ftn58                  = [ftn58;iijjttdd_R13'];

%VR2 = tvary(6)*(C_R21*CxCy+C_R22*CyCz+C_R23*CzCx);
[ii_R21,jj_R21,tt_R21] = find(sparse(C_R21));
tt_R21                 = (1/2)*(1/2)*tvary(6)*tt_R21;
iijjtt_R21             = [ii_R21,jj_R21,tt_R21];
iijjttdd_R21           = combvec(iijjtt_R21',dCxCy');
ftn58                  = [ftn58;iijjttdd_R21'];

[ii_R22,jj_R22,tt_R22] = find(sparse(C_R22));
tt_R22                 = (1/2)*(1/2)*tvary(6)*tt_R22;
iijjtt_R22             = [ii_R22,jj_R22,tt_R22];
iijjttdd_R22           = combvec(iijjtt_R22',dCyCz');
ftn58                  = [ftn58;iijjttdd_R22'];

[ii_R23,jj_R23,tt_R23] = find(sparse(C_R23));
tt_R23                 = (1/2)*(1/2)*tvary(5)*tt_R23;
iijjtt_R23             = [ii_R23,jj_R23,tt_R23];
iijjttdd_R23           = combvec(iijjtt_R23',dCzCx');
ftn58                  = [ftn58;iijjttdd_R23'];

%VR3 = tvary(7)*(C_R31*SxSy+C_R32*SySz+C_R33*SzSx);
dSxSy_p	= [[0 0 0];[-1 -1 0]];
dSxSy_n = [[0 -1 0];[-1 0 0]];
dSySz_p	= [[0 0 0];[0 -1 -1]];
dSySz_n = [[0 -1 0];[0 0 -1]];
dSxSz_p	= [[0 0 0];[-1 0 -1]];
dSxSz_n = [[0 0 -1];[-1 0 0]];

[ii_R31_p,jj_R31_p,tt_R31_p] = find(sparse(C_R31));
tt_R31_p                     = (1/2i)*(1/2i)*tvary(7)*tt_R31_p;
iijjtt_R31_p                 = [ii_R31_p,jj_R31_p,tt_R31_p];
iijjttdd_R31_p     		     = combvec(iijjtt_R31_p',dSxSy_p');
ftn58               	     = [ftn58;iijjttdd_R31_p'];

[ii_R31_n,jj_R31_n,tt_R31_n] = find(sparse(C_R31));
tt_R31_n                     = (-1)*(1/2i)*(1/2i)*tvary(7)*tt_R31_n;
iijjtt_R31_n                 = [ii_R31_n,jj_R31_n,tt_R31_n];
iijjttdd_R31_n     		     = combvec(iijjtt_R31_n',dSxSy_n');
ftn58               	     = [ftn58;iijjttdd_R31_n'];

[ii_R32_p,jj_R32_p,tt_R32_p] = find(sparse(C_R32));
tt_R32_p                     = (1/2i)*(1/2i)*tvary(7)*tt_R32_p;
iijjtt_R32_p                 = [ii_R32_p,jj_R32_p,tt_R32_p];
iijjttdd_R32_p     		     = combvec(iijjtt_R32_p',dSySz_p');
ftn58               	     = [ftn58;iijjttdd_R32_p'];

[ii_R32_n,jj_R32_n,tt_R32_n] = find(sparse(C_R32));
tt_R32_n                     = (-1)*(1/2i)*(1/2i)*tvary(7)*tt_R32_n;
iijjtt_R32_n                 = [ii_R32_n,jj_R32_n,tt_R32_n];
iijjttdd_R32_n     		     = combvec(iijjtt_R32_n',dSySz_n');
ftn58               	     = [ftn58;iijjttdd_R32_n'];

[ii_R33_p,jj_R33_p,tt_R33_p] = find(sparse(C_R33));
tt_R33_p                     = (1/2i)*(1/2i)*tvary(7)*tt_R33_p;
iijjtt_R33_p                 = [ii_R33_p,jj_R33_p,tt_R33_p];
iijjttdd_R33_p     		     = combvec(iijjtt_R33_p',dSxSz_p');
ftn58               	     = [ftn58;iijjttdd_R33_p'];

[ii_R33_n,jj_R33_n,tt_R33_n] = find(sparse(C_R33));
tt_R33_n                     = (-1)*(1/2i)*(1/2i)*tvary(7)*tt_R33_n;
iijjtt_R33_n                 = [ii_R33_n,jj_R33_n,tt_R33_n];
iijjttdd_R33_n     		     = combvec(iijjtt_R33_n',dSxSz_n');
ftn58               	     = [ftn58;iijjttdd_R33_n'];

%VS1 = tvary(8)*(C_S11*SxCy+C_S12*SyCz+C_S13*SzCx);
dSxCy_p = [[0 0 0];[0 -1 0]];
dSxCy_n = [[-1 0 0];[-1 -1 0]];
dSyCz_p = [[0 0 0];[0 0 -1]];
dSyCz_n = [[0 -1 0];[0 -1 -1]];
dSzCx_p = [[0 0 0];[-1 0 0]];
dSzCx_n = [[0 0 -1];[-1 0 -1]];

[ii_S11_p,jj_S11_p,tt_S11_p] = find(sparse(C_S11));
tt_S11_p                     = (1/2i)*(1/2)*tvary(8)*tt_S11_p;
iijjtt_S11_p                 = [ii_S11_p,jj_S11_p,tt_S11_p];
iijjttdd_S11_p     		     = combvec(iijjtt_S11_p',dSxCy_p');
ftn58               	     = [ftn58;iijjttdd_S11_p'];

[ii_S11_n,jj_S11_n,tt_S11_n] = find(sparse(C_S11));
tt_S11_n                     = (-1)*(1/2i)*(1/2)*tvary(7)*tt_S11_n;
iijjtt_S11_n                 = [ii_S11_n,jj_S11_n,tt_S11_n];
iijjttdd_S11_n     		     = combvec(iijjtt_S11_n',dSxCy_n');
ftn58               	     = [ftn58;iijjttdd_S11_n'];

[ii_S12_p,jj_S12_p,tt_S12_p] = find(sparse(C_S12));
tt_S12_p                     = (1/2i)*(1/2)*tvary(8)*tt_S12_p;
iijjtt_S12_p                 = [ii_S12_p,jj_S12_p,tt_S12_p];
iijjttdd_S12_p     		     = combvec(iijjtt_S12_p',dSyCz_p');
ftn58               	     = [ftn58;iijjttdd_S12_p'];

[ii_S12_n,jj_S12_n,tt_S12_n] = find(sparse(C_S12));
tt_S12_n                     = (-1)*(1/2i)*(1/2)*tvary(8)*tt_S12_n;
iijjtt_S12_n                 = [ii_S12_n,jj_S12_n,tt_S12_n];
iijjttdd_S12_n     		     = combvec(iijjtt_S12_n',dSyCz_n');
ftn58               	     = [ftn58;iijjttdd_S12_n'];

[ii_S13_p,jj_S13_p,tt_S13_p] = find(sparse(C_S13));
tt_S13_p                     = (1/2i)*(1/2)*tvary(8)*tt_S13_p;
iijjtt_S13_p                 = [ii_S13_p,jj_S13_p,tt_S13_p];
iijjttdd_S13_p     		     = combvec(iijjtt_S13_p',dSzCx_p');
ftn58               	     = [ftn58;iijjttdd_S13_p'];

[ii_S13_n,jj_S13_n,tt_S13_n] = find(sparse(C_S13));
tt_S13_n                     = (-1)*(1/2i)*(1/2)*tvary(8)*tt_S13_n;
iijjtt_S13_n                 = [ii_S13_n,jj_S13_n,tt_S13_n];
iijjttdd_S13_n     		     = combvec(iijjtt_S13_n',dSzCx_n');
ftn58               	     = [ftn58;iijjttdd_S13_n'];

%VS2 = tvary(9)*(C_S21*CxSy+C_S22*CySz+C_S23*CzSx);
[ii_S21_p,jj_S21_p,tt_S21_p] = find(sparse(C_S21));
tt_S21_p     	             = (1/2)*(1/2i)*tvary(9)*tt_S21_p;
iijjtt_S21_p                 = [ii_S21_p,jj_S21_p,tt_S21_p];
iijjttdd_S21_p     		     = combvec(iijjtt_S21_p',dCxSy_p');
ftn58               	     = [ftn58;iijjttdd_S21_p'];

[ii_S21_n,jj_S21_n,tt_S21_n] = find(sparse(C_S21));
tt_S21_n                     = (-1)*(1/2)*(1/2i)*tvary(9)*tt_S21_n;
iijjtt_S21_n                 = [ii_S21_n,jj_S21_n,tt_S21_n];
iijjttdd_S21_n     		     = combvec(iijjtt_S21_n',dCxSy_n');
ftn58               	     = [ftn58;iijjttdd_S21_n'];

[ii_S22_p,jj_S22_p,tt_S22_p] = find(sparse(C_S22));
tt_S22_p                     = (1/2)*(1/2i)*tvary(9)*tt_S22_p;
iijjtt_S22_p                 = [ii_S22_p,jj_S22_p,tt_S22_p];
iijjttdd_S22_p     		     = combvec(iijjtt_S22_p',dCySz_p');
ftn58               	     = [ftn58;iijjttdd_S22_p'];

[ii_S22_n,jj_S22_n,tt_S22_n] = find(sparse(C_S22));
tt_S22_n                     = (-1)*(1/2)*(1/2i)*tvary(9)*tt_S22_n;
iijjtt_S22_n                 = [ii_S22_n,jj_S22_n,tt_S22_n];
iijjttdd_S22_n               = combvec(iijjtt_S22_n',dCySz_n');
ftn58               	     = [ftn58;iijjttdd_S22_n'];

[ii_S23_p,jj_S23_p,tt_S23_p] = find(sparse(C_S23));
tt_S23_p                     = (1/2)*(1/2i)*tvary(9)*tt_S23_p;
iijjtt_S23_p                 = [ii_S23_p,jj_S23_p,tt_S23_p];
iijjttdd_S23_p     		     = combvec(iijjtt_S23_p',dCzSx_p');
ftn58               	     = [ftn58;iijjttdd_S23_p'];

[ii_S23_n,jj_S23_n,tt_S23_n] = find(sparse(C_S23));
tt_S23_n                     = (-1)*(1/2)*(1/2i)*tvary(9)*tt_S23_n;
iijjtt_S23_n                 = [ii_S23_n,jj_S23_n,tt_S23_n];
iijjttdd_S23_n       	     = combvec(iijjtt_S23_n',dCzSx_n');
ftn58               	     = [ftn58;iijjttdd_S23_n'];

%VS3 = tvary(10)*(C_S31*CxSy+C_S32*CySz+C_S33*CzSx);
[ii_S31_p,jj_S31_p,tt_S31_p] = find(sparse(C_S31));
tt_S31_p     	             = (1/2)*(1/2i)*tvary(10)*tt_S31_p;
iijjtt_S31_p                 = [ii_S31_p,jj_S31_p,tt_S31_p];
iijjttdd_S31_p     		     = combvec(iijjtt_S31_p',dCxSy_p');
ftn58               	     = [ftn58;iijjttdd_S31_p'];

[ii_S31_n,jj_S31_n,tt_S31_n] = find(sparse(C_S31));
tt_S31_n                     = (-1)*(1/2)*(1/2i)*tvary(10)*tt_S31_n;
iijjtt_S31_n                 = [ii_S31_n,jj_S31_n,tt_S31_n];
iijjttdd_S31_n     		     = combvec(iijjtt_S31_n',dCxSy_n');
ftn58               	     = [ftn58;iijjttdd_S31_n'];

[ii_S32_p,jj_S32_p,tt_S32_p] = find(sparse(C_S32));
tt_S32_p                     = (1/2)*(1/2i)*tvary(10)*tt_S32_p;
iijjtt_S32_p                 = [ii_S32_p,jj_S32_p,tt_S32_p];
iijjttdd_S32_p     		     = combvec(iijjtt_S32_p',dCySz_p');
ftn58               	     = [ftn58;iijjttdd_S32_p'];

[ii_S32_n,jj_S32_n,tt_S32_n] = find(sparse(C_S32));
tt_S32_n                     = (-1)*(1/2)*(1/2i)*tvary(10)*tt_S32_n;
iijjtt_S32_n                 = [ii_S32_n,jj_S32_n,tt_S32_n];
iijjttdd_S32_n               = combvec(iijjtt_S32_n',dCySz_n');
ftn58               	     = [ftn58;iijjttdd_S32_n'];

[ii_S33_p,jj_S33_p,tt_S33_p] = find(sparse(C_S33));
tt_S33_p                     = (1/2)*(1/2i)*tvary(10)*tt_S33_p;
iijjtt_S33_p                 = [ii_S33_p,jj_S33_p,tt_S33_p];
iijjttdd_S33_p     		     = combvec(iijjtt_S33_p',dCzSx_p');
ftn58               	     = [ftn58;iijjttdd_S33_p'];

[ii_S33_n,jj_S33_n,tt_S33_n] = find(sparse(C_S33));
tt_S33_n                     = (-1)*(1/2)*(1/2i)*tvary(10)*tt_S33_n;
iijjtt_S33_n                 = [ii_S33_n,jj_S33_n,tt_S33_n];
iijjttdd_S33_n       	     = combvec(iijjtt_S33_n',dCzSx_n');
ftn58               	     = [ftn58;iijjttdd_S33_n'];

%==========================================================================
%ij = ftn58(2:end,1:2);
%tt = ftn58(2:end,3);
%dd = ftn58(2:end,4:6);

%% compensate the lower triangular part
%v          = find(ij(:,1)==ij(:,2));
%ij_temp    = ij;
%tt_temp    = tt;
%dd_temp    = dd;
%ij_temp(v,:) = [];
%tt_temp(v) = []; 
%dd_temp(v,:) = [];
%ftn58      = [ftn58;[flip(ij_temp,2) conj(tt_temp) (-1)*dd_temp]];

ftn58sparse.ij = ftn58(2:end,1:2);
ftn58sparse.tt = ftn58(2:end,3);
ftn58sparse.dd = ftn58(2:end,4:6);

save RhSi_WSOC_ftn58sparse_ver6.mat ftn58sparse

end