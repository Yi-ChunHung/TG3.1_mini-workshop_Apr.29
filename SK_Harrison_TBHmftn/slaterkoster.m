%%% Slater-Koster angle denpendence Tight-Binding model
%%%---------------------------------------------------------------------------------------------------------------
%%% origin model : J. C. SLATER AND G. F. KOSTER, Phys. Rev. 94, 1498 (1954), Table 1
%%% easier article : W. E. Pickett, “Tight Binding” Method: Linear Combination of Atomic Orbitals (LCAO) (2014)
%%
function [ESK] = slaterkoster(l,m,n,ob1,ob2)
% l,m,n: direction cosines (for x,y,z axis,respectively)
% 1, 2, 3, 4, 5,  6,  7,  8,    9
% s,px,py,pz,dz,dxz,dyz,dxy,dx2y2
global sss sps pps ppp sds pds pdp dds ddp ddd

%% make sure that only the upper triangular elements are constructed
if ob1>ob2
    obtemp = ob2;
    ob2    = ob1;
    ob1    = obtemp;
    l=-l;    m=-m;    n=-n;
end

%% generate hopping strength through SK method
if ob1==1 & ob2==1
    ESK=sss;
    return
elseif ob1==1 & ob2==2
    ESK=espx(l);
    return
elseif ob1==1 & ob2==3
    ESK=espx(m);
    return
elseif ob1==1 & ob2==4
    ESK=espx(n);
    return
elseif ob1==2 & ob2==2
    ESK=epxpx(l);
    return
elseif ob1==3 & ob2==3
    ESK=epxpx(m);
    return
elseif ob1==4 & ob2==4
    ESK=epxpx(n);
    return
elseif ob1==2 & ob2==3
    ESK=epxpy(l,m);
    return
elseif ob1==2 & ob2==4
    ESK=epxpy(l,n);
    return
elseif ob1==3 & ob2==4
    ESK=epxpy(m,n);
    return
elseif ob1==1 & ob2==5
    ESK=(n*n-(l*l+m*m)/2)*sds;
    return
elseif ob1==1 & ob2==6
    ESK=esdxy(l,n);
    return
elseif ob1==1 & ob2==7
    ESK=esdxy(m,n);
    return
elseif ob1==1 & ob2==8
    ESK=esdxy(l,m);
    return
elseif ob1==1 & ob2==9
    ESK=sqrt(3)/2*(l*l-m*m)*sds;
    return
elseif ob1==2 & ob2==5
    ESK=epxdz(l,m,n);
    return
elseif ob1==3 & ob2==5
    ESK=epxdz(m,l,n);
    return
elseif ob1==4 & ob2==5
    ESK=n*(n*n-(l*l+m*m)/2)*pds+sqrt(3)*n*(l*l+m*m)*pdp;
    return
elseif ob1==2 & ob2==6
    ESK=epxdxy(l,n);
    return
elseif ob1==2 & ob2==8
    ESK=epxdxy(l,m);
    return
elseif ob1==3 & ob2==7
    ESK=epxdxy(m,n);
    return
elseif ob1==3 & ob2==8
    ESK=epxdxy(m,l);
    return
elseif ob1==4 & ob2==6
    ESK=epxdxy(n,l);
    return
elseif ob1==4 & ob2==7
    ESK=epxdxy(n,m);
    return
elseif ob1==2 & ob2==7
    ESK=l*m*n*(sqrt(3)*pds-2*pdp);
    return
elseif ob1==3 & ob2==6
    ESK=l*m*n*(sqrt(3)*pds-2*pdp);
    return
elseif ob1==4 & ob2==8
    ESK=l*m*n*(sqrt(3)*pds-2*pdp);
    return
elseif ob1==2 & ob2==9
    ESK=epxdx2y2(l,m,n);
    return
elseif ob1==3 & ob2==9
    ESK=-epxdx2y2(m,l,n);
    return
elseif ob1==4 & ob2==9
    ESK=n*(l*l-m*m)*(sqrt(3)/2*pds-pdp);
    return
elseif ob1==5 & ob2==5
    ESK=(n*n-(l*l+m*m)/2)^2*dds+3*n*n*(l*l+m*m)*ddp+3/4*(l*l+m*m)^2*ddd;
    return
elseif ob1==5 & ob2==6
    ESK=edzdxz(l,m,n);
    return
elseif ob1==5 & ob2==7
    ESK=edzdxz(m,l,n);
    return
elseif ob1==5 & ob2==8
    ESK=sqrt(3)*l*m*((n*n-(l*l+m*m)/2)*dds-2*n*n*ddp+(1+n*n)/2*ddd);
    return
elseif ob1==5 & ob2==9
    ESK=sqrt(3)*(l*l-m*m)*((n*n-(l*l+m*m)/2)/2*dds-n*n*ddp+(1+n*n)/4*ddd);
    return
elseif ob1==6 & ob2==6
    ESK=edxydxy(l,n,m);
    return
elseif ob1==7 & ob2==7
    ESK=edxydxy(m,n,l);
    return
elseif ob1==8 & ob2==8
    ESK=edxydxy(l,m,n);
    return
elseif ob1==6 & ob2==7 %dxz dyz
    ESK=edxydxz(n,l,m);
    return
elseif ob1==6 & ob2==8 %dxz dxy
    ESK=edxydxz(l,m,n);
    return
elseif ob1==7 & ob2==8 %dyz dxy
    ESK=edxydxz(m,n,l);
    return
elseif ob1==6 & ob2==9 %dxz dx2y2
    ESK=edxzdx2y2(l,m,n);
    return
elseif ob1==7 & ob2==9 %dyz dx2y2
    ESK=-edxzdx2y2(m,l,n);
    return
elseif ob1==8 & ob2==9 %dxy dx2y2
    ESK=l*m*(3/2*(l*l-m*m)*dds+2*(m*m-l*l)*ddp+(l*l-m*m)/2*ddd);
    return
elseif ob1==9 & ob2==9 %dx2y2
    ESK=3/4*(l*l-m*m)^2*dds+(l*l+m*m-(l*l-m*m)^2)*ddp+(n*n+(l*l-m*m)^2/4)*ddd;
    return   
else
    error('slaterkoster:orbital wrong');
end

%% -------------------------------------------------------------------- %%
%% SK methods
function [esp]=espx(ll)
global sps
esp=ll*sps;
return

function [esp]=epxpx(ll)
global pps ppp
esp=ll*ll*pps+(1-ll*ll)*ppp;
return

function [esp]=epxpy(ll,mm)
global pps ppp
esp=ll*mm*(pps-ppp);
return

function [esp]=esdxy(ll,mm)
global sds
esp=sqrt(3)*ll*mm*sds;
return

function [esp]=epxdz(ll,mm,nn)
global pds pdp
esp=ll*(nn*nn-(ll*ll+mm*mm)/2)*pds-sqrt(3)*ll*nn*nn*pdp;
return

function [esp]=epxdx2y2(ll,mm,nn)
global pds pdp
esp=sqrt(3)/2*ll*(ll*ll-mm*mm)*pds+ll*(1-ll*ll+mm*mm)*pdp;
return

function [esp]=epxdxy(ll,mm)
global pds pdp
esp=sqrt(3)*ll*ll*mm*pds+mm*(1-2*ll*ll)*pdp;
return

function [esp]=edzdxz(ll,mm,nn)
global dds ddp ddd
esp=sqrt(3)*ll*nn*((nn*nn-(ll*ll+mm*mm)/2)*dds+(ll*ll+mm*mm-nn*nn)*ddp-(ll*ll+mm*mm)/2*ddd);
return

function [esp]=edxydxy(ll,mm,nn)
global dds ddp ddd
esp=3*ll*ll*mm*mm*dds+(ll*ll+mm*mm-4*ll*ll*mm*mm)*ddp+(nn*nn+ll*ll*mm*mm)*ddd;
return

function [esp]=edxydxz(ll,mm,nn)
global dds ddp ddd
esp=mm*nn*(3*ll*ll*dds+(1-4*ll*ll)*ddp+(ll*ll-1)*ddd);
return

function [esp]=edxzdx2y2(ll,mm,nn)
global dds ddp ddd
esp=nn*ll*(3/2*(ll*ll-mm*mm)*dds+(1-2*(ll*ll-mm*mm))*ddp-(1-(ll*ll-mm*mm)/2)*ddd);
return
