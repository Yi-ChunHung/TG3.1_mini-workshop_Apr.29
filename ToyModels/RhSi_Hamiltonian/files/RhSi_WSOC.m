classdef RhSi_WSOC 
    methods(Static)
%==========================================================================
%                   --  RhSi Hamiltonian With  SOC --              
%==========================================================================

%==========================================================================
%                               -- H --
%==========================================================================
%==========================================================================
function H = RhSi_H(kp)
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
CxCy = cos(kp(1,1)/2)*cos(kp(1,2)/2);
CyCz = cos(kp(1,2)/2)*cos(kp(1,3)/2);
CzCx = cos(kp(1,1)/2)*cos(kp(1,3)/2);
    
CxSy = cos(kp(1,1)/2)*sin(kp(1,2)/2);
CySz = cos(kp(1,2)/2)*sin(kp(1,3)/2);
CzSx = cos(kp(1,3)/2)*sin(kp(1,1)/2);
    
SxSy = sin(kp(1,1)/2)*sin(kp(1,2)/2);
SySz = sin(kp(1,2)/2)*sin(kp(1,3)/2);
SzSx = sin(kp(1,1)/2)*sin(kp(1,3)/2);
    
SxCy = sin(kp(1,1)/2)*cos(kp(1,2)/2);
SyCz = sin(kp(1,2)/2)*cos(kp(1,3)/2);
SzCx = sin(kp(1,3)/2)*cos(kp(1,1)/2);
%--------------------------------------------------------------------------
 
%--------------------------------------------------------------------------
H1  = tvary(1)*(C11*CxCy+C12*CyCz+C13*CzCx);
H2  = tvary(2)*(C21*CxSy+C22*CySz+C23*CzSx);

H3  = tvary(3)*C3*(cos(kp(1,1))+cos(kp(1,2))+cos(kp(1,3)));

VR1 = tvary(5)*(C_R11*CxCy+C_R12*CyCz+C_R13*CzCx);
VR2 = tvary(6)*(C_R21*CxCy+C_R22*CyCz+C_R23*CzCx);
VR3 = tvary(7)*(C_R31*SxSy+C_R32*SySz+C_R33*SzSx);
VS1 = tvary(8)*(C_S11*SxCy+C_S12*SyCz+C_S13*SzCx);
VS2 = tvary(9)*(C_S21*CxSy+C_S22*CySz+C_S23*CzSx);
VS3 = tvary(10)*(C_S31*CxSy+C_S32*CySz+C_S33*CzSx);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
H = H1+H2+H3+VR1+VR2+VR3+VS1+VS2+VS3; 
%--------------------------------------------------------------------------

end
%==========================================================================

%==========================================================================
%                              --  dH --              
%==========================================================================
%==========================================================================
function dHx = RhSi_dHx(kp)
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
CxCy = -1/2*sin(kp(1,1)/2)   *       cos(kp(1,2)/2);
CyCz = 0;
CzCx = -1/2*sin(kp(1,1)/2)   *       cos(kp(1,3)/2);
    
CxSy = -1/2*sin(kp(1,1)/2)   *       sin(kp(1,2)/2);
CySz = 0;
CzSx =      cos(kp(1,3)/2)   *   1/2*cos(kp(1,1)/2);
    
SxSy = 1/2*cos(kp(1,1)/2)    *       sin(kp(1,2)/2);
SySz = 0;
SzSx = 1/2*cos(kp(1,1)/2)    *       sin(kp(1,3)/2);
    
SxCy = 1/2*cos(kp(1,1)/2)    *        cos(kp(1,2)/2);
SyCz = 0;
SzCx =     sin(kp(1,3)/2)    *   -1/2*sin(kp(1,1)/2);
%--------------------------------------------------------------------------
 
%--------------------------------------------------------------------------
H1  = tvary(1)*(C11*CxCy+C12*CyCz+C13*CzCx);
H2  = tvary(2)*(C21*CxSy+C22*CySz+C23*CzSx);

H3  = tvary(3)*C3*(-sin(kp(1,1)));

VR1 = tvary(5)*(C_R11*CxCy+C_R12*CyCz+C_R13*CzCx);
VR2 = tvary(6)*(C_R21*CxCy+C_R22*CyCz+C_R23*CzCx);
VR3 = tvary(7)*(C_R31*SxSy+C_R32*SySz+C_R33*SzSx);
VS1 = tvary(8)*(C_S11*SxCy+C_S12*SyCz+C_S13*SzCx);
VS2 = tvary(9)*(C_S21*CxSy+C_S22*CySz+C_S23*CzSx);
VS3 = tvary(10)*(C_S31*CxSy+C_S32*CySz+C_S33*CzSx);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
dHx = H1+H2+H3+VR1+VR2+VR3+VS1+VS2+VS3; 
%--------------------------------------------------------------------------

end
%==========================================================================

%==========================================================================
function dHy = RhSi_dHy(kp)
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
CxCy =      cos(kp(1,1)/2)   *  -1/2*sin(kp(1,2)/2);
CyCz = -1/2*sin(kp(1,2)/2)   *       cos(kp(1,3)/2);
CzCx = 0;
    
CxSy =      cos(kp(1,1)/2)   *   1/2*cos(kp(1,2)/2);
CySz = -1/2*sin(kp(1,2)/2)   *       sin(kp(1,3)/2);
CzSx = 0;
    
SxSy =      sin(kp(1,1)/2)   *   1/2*cos(kp(1,2)/2);
SySz =  1/2*cos(kp(1,2)/2)   *       sin(kp(1,3)/2);
SzSx = 0;
    
SxCy =      sin(kp(1,1)/2)   *  -1/2*sin(kp(1,2)/2);
SyCz =  1/2*cos(kp(1,2)/2)   *       cos(kp(1,3)/2);
SzCx = 0;
%--------------------------------------------------------------------------
 
%--------------------------------------------------------------------------
H1  = tvary(1)*(C11*CxCy+C12*CyCz+C13*CzCx);
H2  = tvary(2)*(C21*CxSy+C22*CySz+C23*CzSx);

H3  = tvary(3)*C3*(-sin(kp(1,2)));

VR1 = tvary(5)*(C_R11*CxCy+C_R12*CyCz+C_R13*CzCx);
VR2 = tvary(6)*(C_R21*CxCy+C_R22*CyCz+C_R23*CzCx);
VR3 = tvary(7)*(C_R31*SxSy+C_R32*SySz+C_R33*SzSx);
VS1 = tvary(8)*(C_S11*SxCy+C_S12*SyCz+C_S13*SzCx);
VS2 = tvary(9)*(C_S21*CxSy+C_S22*CySz+C_S23*CzSx);
VS3 = tvary(10)*(C_S31*CxSy+C_S32*CySz+C_S33*CzSx);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
dHy = H1+H2+H3+VR1+VR2+VR3+VS1+VS2+VS3; 
%--------------------------------------------------------------------------
end
%==========================================================================

%==========================================================================
function dHz = RhSi_dHz(kp)
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
CxCy = 0;
CyCz =      cos(kp(1,2)/2)      * -1/2*sin(kp(1,3)/2);
CzCx =      cos(kp(1,1)/2)      * -1/2*sin(kp(1,3)/2);
    
CxSy = 0;
CySz =      cos(kp(1,2)/2)      *  1/2*cos(kp(1,3)/2);
CzSx = -1/2*sin(kp(1,3)/2)      *      sin(kp(1,1)/2);
    
SxSy = 0;
SySz =      sin(kp(1,2)/2)      *  1/2*cos(kp(1,3)/2);
SzSx =      sin(kp(1,1)/2)      *  1/2*cos(kp(1,3)/2);
    
SxCy = 0;
SyCz =      sin(kp(1,2)/2)      * -1/2*sin(kp(1,3)/2);
SzCx =  1/2*cos(kp(1,3)/2)      *      cos(kp(1,1)/2);
%--------------------------------------------------------------------------
 
%--------------------------------------------------------------------------
H1  = tvary(1)*(C11*CxCy+C12*CyCz+C13*CzCx);
H2  = tvary(2)*(C21*CxSy+C22*CySz+C23*CzSx);

H3  = tvary(3)*C3*(-sin(kp(1,3)));

VR1 = tvary(5)*(C_R11*CxCy+C_R12*CyCz+C_R13*CzCx);
VR2 = tvary(6)*(C_R21*CxCy+C_R22*CyCz+C_R23*CzCx);
VR3 = tvary(7)*(C_R31*SxSy+C_R32*SySz+C_R33*SzSx);
VS1 = tvary(8)*(C_S11*SxCy+C_S12*SyCz+C_S13*SzCx);
VS2 = tvary(9)*(C_S21*CxSy+C_S22*CySz+C_S23*CzSx);
VS3 = tvary(10)*(C_S31*CxSy+C_S32*CySz+C_S33*CzSx);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
dHz = H1+H2+H3+VR1+VR2+VR3+VS1+VS2+VS3; 
%--------------------------------------------------------------------------
end

%==========================================================================


%==========================================================================
%                              --  dH dH --              
%==========================================================================
%==========================================================================
function dHx_dHx = RhSi_dHx_dHx(kp)
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
CxCy = -1/4*cos(kp(1,1)/2)   *       cos(kp(1,2)/2);
CyCz = 0;
CzCx = -1/4*cos(kp(1,1)/2)   *       cos(kp(1,3)/2);
    
CxSy = -1/4*cos(kp(1,1)/2)   *       sin(kp(1,2)/2);
CySz = 0;
CzSx =      cos(kp(1,3)/2)   *  -1/4*sin(kp(1,1)/2);
    
SxSy = -1/4*sin(kp(1,1)/2)   *       sin(kp(1,2)/2);
SySz = 0;
SzSx = -1/4*sin(kp(1,1)/2)   *       sin(kp(1,3)/2);
    
SxCy = -1/4*sin(kp(1,1)/2)   *       cos(kp(1,2)/2);
SyCz = 0;
SzCx =      sin(kp(1,3)/2)   *  -1/4*cos(kp(1,1)/2);
%--------------------------------------------------------------------------
 
%--------------------------------------------------------------------------
H1  = tvary(1)*(C11*CxCy+C12*CyCz+C13*CzCx);
H2  = tvary(2)*(C21*CxSy+C22*CySz+C23*CzSx);

H3  = tvary(3)*C3*(-cos(kp(1,1)));

VR1 = tvary(5)*(C_R11*CxCy+C_R12*CyCz+C_R13*CzCx);
VR2 = tvary(6)*(C_R21*CxCy+C_R22*CyCz+C_R23*CzCx);
VR3 = tvary(7)*(C_R31*SxSy+C_R32*SySz+C_R33*SzSx);
VS1 = tvary(8)*(C_S11*SxCy+C_S12*SyCz+C_S13*SzCx);
VS2 = tvary(9)*(C_S21*CxSy+C_S22*CySz+C_S23*CzSx);
VS3 = tvary(10)*(C_S31*CxSy+C_S32*CySz+C_S33*CzSx);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
dHx_dHx = H1+H2+H3+VR1+VR2+VR3+VS1+VS2+VS3; 
%--------------------------------------------------------------------------

end
%==========================================================================


%==========================================================================
function dHx_dHy = RhSi_dHx_dHy(kp)
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
CxCy = -1/2*sin(kp(1,1)/2)   *  -1/2*sin(kp(1,2)/2);
CyCz = 0;
CzCx = 0;
    
CxSy = -1/2*sin(kp(1,1)/2)   *   1/2*cos(kp(1,2)/2);
CySz = 0;
CzSx = 0;
    
SxSy = 1/2*cos(kp(1,1)/2)    *   1/2*cos(kp(1,2)/2);
SySz = 0;
SzSx = 0;
    
SxCy = 1/2*cos(kp(1,1)/2)    *   -1/2*sin(kp(1,2)/2);
SyCz = 0;
SzCx = 0;
 %--------------------------------------------------------------------------
 
 %--------------------------------------------------------------------------
H1  = tvary(1)*(C11*CxCy+C12*CyCz+C13*CzCx);
H2  = tvary(2)*(C21*CxSy+C22*CySz+C23*CzSx);

H3  = tvary(3)*C3*(0);

VR1 = tvary(5)*(C_R11*CxCy+C_R12*CyCz+C_R13*CzCx);
VR2 = tvary(6)*(C_R21*CxCy+C_R22*CyCz+C_R23*CzCx);
VR3 = tvary(7)*(C_R31*SxSy+C_R32*SySz+C_R33*SzSx);
VS1 = tvary(8)*(C_S11*SxCy+C_S12*SyCz+C_S13*SzCx);
VS2 = tvary(9)*(C_S21*CxSy+C_S22*CySz+C_S23*CzSx);
VS3 = tvary(10)*(C_S31*CxSy+C_S32*CySz+C_S33*CzSx);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
dHx_dHy = H1+H2+H3+VR1+VR2+VR3+VS1+VS2+VS3; 
%--------------------------------------------------------------------------

end
%==========================================================================


%==========================================================================
function dHx_dHz = RhSi_dHx_dHz(kp)
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
CxCy = 0;
CyCz = 0;
CzCx = -1/2*sin(kp(1,1)/2)   *  -1/2*sin(kp(1,3)/2);
    
CxSy = 0;
CySz = 0;
CzSx = -1/2*sin(kp(1,3)/2)   *   1/2*cos(kp(1,1)/2);
    
SxSy = 0;
SySz = 0;
SzSx = 1/2*cos(kp(1,1)/2)    *   1/2*cos(kp(1,3)/2);
    
SxCy = 0;
SyCz = 0;
SzCx = 1/2*cos(kp(1,3)/2)    *   -1/2*sin(kp(1,1)/2);
 %--------------------------------------------------------------------------
 
 %--------------------------------------------------------------------------
H1  = tvary(1)*(C11*CxCy+C12*CyCz+C13*CzCx);
H2  = tvary(2)*(C21*CxSy+C22*CySz+C23*CzSx);

H3  = tvary(3)*C3*(0);

VR1 = tvary(5)*(C_R11*CxCy+C_R12*CyCz+C_R13*CzCx);
VR2 = tvary(6)*(C_R21*CxCy+C_R22*CyCz+C_R23*CzCx);
VR3 = tvary(7)*(C_R31*SxSy+C_R32*SySz+C_R33*SzSx);
VS1 = tvary(8)*(C_S11*SxCy+C_S12*SyCz+C_S13*SzCx);
VS2 = tvary(9)*(C_S21*CxSy+C_S22*CySz+C_S23*CzSx);
VS3 = tvary(10)*(C_S31*CxSy+C_S32*CySz+C_S33*CzSx);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
dHx_dHz = H1+H2+H3+VR1+VR2+VR3+VS1+VS2+VS3; 
%--------------------------------------------------------------------------

end
%==========================================================================


%==========================================================================
function dHy_dHx = RhSi_dHy_dHx(kp)
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
CxCy = -1/2*sin(kp(1,1)/2)   *  -1/2*sin(kp(1,2)/2);
CyCz = 0;
CzCx = 0;
    
CxSy = -1/2*sin(kp(1,1)/2)   *   1/2*cos(kp(1,2)/2);
CySz = 0;
CzSx = 0;
    
SxSy = 1/2*cos(kp(1,1)/2)    *   1/2*cos(kp(1,2)/2);
SySz = 0;
SzSx = 0;
    
SxCy = 1/2*cos(kp(1,1)/2)    *   -1/2*sin(kp(1,2)/2);
SyCz = 0;
SzCx = 0;
 %--------------------------------------------------------------------------
 
 %--------------------------------------------------------------------------
H1  = tvary(1)*(C11*CxCy+C12*CyCz+C13*CzCx);
H2  = tvary(2)*(C21*CxSy+C22*CySz+C23*CzSx);

H3  = tvary(3)*C3*(0);

VR1 = tvary(5)*(C_R11*CxCy+C_R12*CyCz+C_R13*CzCx);
VR2 = tvary(6)*(C_R21*CxCy+C_R22*CyCz+C_R23*CzCx);
VR3 = tvary(7)*(C_R31*SxSy+C_R32*SySz+C_R33*SzSx);
VS1 = tvary(8)*(C_S11*SxCy+C_S12*SyCz+C_S13*SzCx);
VS2 = tvary(9)*(C_S21*CxSy+C_S22*CySz+C_S23*CzSx);
VS3 = tvary(10)*(C_S31*CxSy+C_S32*CySz+C_S33*CzSx);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
dHy_dHx = H1+H2+H3+VR1+VR2+VR3+VS1+VS2+VS3; 
%--------------------------------------------------------------------------

end
%==========================================================================


%==========================================================================
function dHy_dHy = RhSi_dHy_dHy(kp)
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
CxCy =      cos(kp(1,1)/2)   *  -1/4*cos(kp(1,2)/2);
CyCz = -1/4*cos(kp(1,2)/2)   *       cos(kp(1,3)/2);
CzCx = 0;
    
CxSy =      cos(kp(1,1)/2)   *  -1/4*sin(kp(1,2)/2);
CySz = -1/4*cos(kp(1,2)/2)   *       sin(kp(1,3)/2);
CzSx = 0;
    
SxSy =      sin(kp(1,1)/2)   *  -1/4*sin(kp(1,2)/2);
SySz = -1/4*sin(kp(1,2)/2)   *       sin(kp(1,3)/2);
SzSx = 0;
    
SxCy =      sin(kp(1,1)/2)   *  -1/4*cos(kp(1,2)/2);
SyCz = -1/4*sin(kp(1,2)/2)   *       cos(kp(1,3)/2);
SzCx = 0;
 %--------------------------------------------------------------------------
 
 %--------------------------------------------------------------------------
H1  = tvary(1)*(C11*CxCy+C12*CyCz+C13*CzCx);
H2  = tvary(2)*(C21*CxSy+C22*CySz+C23*CzSx);

H3  = tvary(3)*C3*(-cos(kp(1,2)));

VR1 = tvary(5)*(C_R11*CxCy+C_R12*CyCz+C_R13*CzCx);
VR2 = tvary(6)*(C_R21*CxCy+C_R22*CyCz+C_R23*CzCx);
VR3 = tvary(7)*(C_R31*SxSy+C_R32*SySz+C_R33*SzSx);
VS1 = tvary(8)*(C_S11*SxCy+C_S12*SyCz+C_S13*SzCx);
VS2 = tvary(9)*(C_S21*CxSy+C_S22*CySz+C_S23*CzSx);
VS3 = tvary(10)*(C_S31*CxSy+C_S32*CySz+C_S33*CzSx);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
dHy_dHy = H1+H2+H3+VR1+VR2+VR3+VS1+VS2+VS3; 
%--------------------------------------------------------------------------

end

%==========================================================================


%==========================================================================
function dHy_dHz = RhSi_dHy_dHz(kp)
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
CxCy = 0;
CyCz = -1/2*sin(kp(1,2)/2)   *  -1/2*sin(kp(1,3)/2);
CzCx = 0;
    
CxSy = 0;
CySz = -1/2*sin(kp(1,2)/2)   *   1/2*cos(kp(1,3)/2);
CzSx = 0;
    
SxSy = 0;
SySz =  1/2*cos(kp(1,2)/2)   *   1/2*cos(kp(1,3)/2);
SzSx = 0;
    
SxCy =  0;
SyCz =  1/2*cos(kp(1,2)/2)   *  -1/2*sin(kp(1,3)/2);
SzCx =  0;
 %--------------------------------------------------------------------------
 
 %--------------------------------------------------------------------------
H1  = tvary(1)*(C11*CxCy+C12*CyCz+C13*CzCx);
H2  = tvary(2)*(C21*CxSy+C22*CySz+C23*CzSx);

H3  = tvary(3)*C3*(0);

VR1 = tvary(5)*(C_R11*CxCy+C_R12*CyCz+C_R13*CzCx);
VR2 = tvary(6)*(C_R21*CxCy+C_R22*CyCz+C_R23*CzCx);
VR3 = tvary(7)*(C_R31*SxSy+C_R32*SySz+C_R33*SzSx);
VS1 = tvary(8)*(C_S11*SxCy+C_S12*SyCz+C_S13*SzCx);
VS2 = tvary(9)*(C_S21*CxSy+C_S22*CySz+C_S23*CzSx);
VS3 = tvary(10)*(C_S31*CxSy+C_S32*CySz+C_S33*CzSx);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
dHy_dHz = H1+H2+H3+VR1+VR2+VR3+VS1+VS2+VS3; 
%--------------------------------------------------------------------------

end

%==========================================================================


%==========================================================================
function dHz_dHx = RhSi_dHz_dHx(kp)
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
CxCy = 0;
CyCz = 0;
CzCx = -1/2*sin(kp(1,1)/2)   *  -1/2*sin(kp(1,3)/2);
    
CxSy = 0;
CySz = 0;
CzSx = -1/2*sin(kp(1,3)/2)   *   1/2*cos(kp(1,1)/2);
    
SxSy = 0;
SySz = 0;
SzSx = 1/2*cos(kp(1,1)/2)    *   1/2*cos(kp(1,3)/2);
    
SxCy = 0;
SyCz = 0;
SzCx = 1/2*cos(kp(1,3)/2)    *   -1/2*sin(kp(1,1)/2);
 %--------------------------------------------------------------------------
 
 %--------------------------------------------------------------------------
H1  = tvary(1)*(C11*CxCy+C12*CyCz+C13*CzCx);
H2  = tvary(2)*(C21*CxSy+C22*CySz+C23*CzSx);

H3  = tvary(3)*C3*(0);

VR1 = tvary(5)*(C_R11*CxCy+C_R12*CyCz+C_R13*CzCx);
VR2 = tvary(6)*(C_R21*CxCy+C_R22*CyCz+C_R23*CzCx);
VR3 = tvary(7)*(C_R31*SxSy+C_R32*SySz+C_R33*SzSx);
VS1 = tvary(8)*(C_S11*SxCy+C_S12*SyCz+C_S13*SzCx);
VS2 = tvary(9)*(C_S21*CxSy+C_S22*CySz+C_S23*CzSx);
VS3 = tvary(10)*(C_S31*CxSy+C_S32*CySz+C_S33*CzSx);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
dHz_dHx = H1+H2+H3+VR1+VR2+VR3+VS1+VS2+VS3; 
%--------------------------------------------------------------------------

end

%==========================================================================


%==========================================================================
function dHz_dHy = RhSi_dHz_dHy(kp)
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
CxCy = 0;
CyCz = -1/2*sin(kp(1,2)/2)   *  -1/2*sin(kp(1,3)/2);
CzCx = 0;
    
CxSy = 0;
CySz = -1/2*sin(kp(1,2)/2)   *   1/2*cos(kp(1,3)/2);
CzSx = 0;
    
SxSy = 0;
SySz =  1/2*cos(kp(1,2)/2)   *   1/2*cos(kp(1,3)/2);
SzSx = 0;
    
SxCy = 0;
SyCz =  1/2*cos(kp(1,2)/2)   *  -1/2*sin(kp(1,3)/2);
SzCx = 0;
 %--------------------------------------------------------------------------
 
 %--------------------------------------------------------------------------
H1  = tvary(1)*(C11*CxCy+C12*CyCz+C13*CzCx);
H2  = tvary(2)*(C21*CxSy+C22*CySz+C23*CzSx);

H3  = tvary(3)*C3*(0);

VR1 = tvary(5)*(C_R11*CxCy+C_R12*CyCz+C_R13*CzCx);
VR2 = tvary(6)*(C_R21*CxCy+C_R22*CyCz+C_R23*CzCx);
VR3 = tvary(7)*(C_R31*SxSy+C_R32*SySz+C_R33*SzSx);
VS1 = tvary(8)*(C_S11*SxCy+C_S12*SyCz+C_S13*SzCx);
VS2 = tvary(9)*(C_S21*CxSy+C_S22*CySz+C_S23*CzSx);
VS3 = tvary(10)*(C_S31*CxSy+C_S32*CySz+C_S33*CzSx);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
dHz_dHy = H1+H2+H3+VR1+VR2+VR3+VS1+VS2+VS3; 
%--------------------------------------------------------------------------

end

%==========================================================================

%==========================================================================
function dHz_dHz = RhSi_dHz_dHz(kp)
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
CxCy = 0;
CyCz =      cos(kp(1,2)/2)      * -1/4*cos(kp(1,3)/2);
CzCx =      cos(kp(1,1)/2)      * -1/4*cos(kp(1,3)/2);
    
CxSy = 0;
CySz =      cos(kp(1,2)/2)      * -1/4*sin(kp(1,3)/2);
CzSx = -1/4*cos(kp(1,3)/2)      *      sin(kp(1,1)/2);
    
SxSy = 0;
SySz =      sin(kp(1,2)/2)      * -1/4*sin(kp(1,3)/2);
SzSx =      sin(kp(1,1)/2)      * -1/4*sin(kp(1,3)/2);
    
SxCy = 0;
SyCz =      sin(kp(1,2)/2)      * -1/4*cos(kp(1,3)/2);
SzCx = -1/4*sin(kp(1,3)/2)      *      cos(kp(1,1)/2);
%--------------------------------------------------------------------------
 
%--------------------------------------------------------------------------
H1  = tvary(1)*(C11*CxCy+C12*CyCz+C13*CzCx);
H2  = tvary(2)*(C21*CxSy+C22*CySz+C23*CzSx);

H3  = tvary(3)*C3*(-cos(kp(1,3)));

VR1 = tvary(5)*(C_R11*CxCy+C_R12*CyCz+C_R13*CzCx);
VR2 = tvary(6)*(C_R21*CxCy+C_R22*CyCz+C_R23*CzCx);
VR3 = tvary(7)*(C_R31*SxSy+C_R32*SySz+C_R33*SzSx);
VS1 = tvary(8)*(C_S11*SxCy+C_S12*SyCz+C_S13*SzCx);
VS2 = tvary(9)*(C_S21*CxSy+C_S22*CySz+C_S23*CzSx);
VS3 = tvary(10)*(C_S31*CxSy+C_S32*CySz+C_S33*CzSx);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
dHz_dHz = H1+H2+H3+VR1+VR2+VR3+VS1+VS2+VS3; 
%--------------------------------------------------------------------------

end
%==========================================================================




    end
end