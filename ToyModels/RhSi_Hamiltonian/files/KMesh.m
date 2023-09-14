
function [kpoints,nks,list] = KMesh(HSP,nk)

SI = size(HSP);

kpoints = [];
list    = [1];
 
for ss = 1:SI(1)-1
 p=[];
 
 if ss ~=SI(1)-1
     
  p = [linspace(HSP(ss,1),HSP(ss+1,1),nk+1)' linspace(HSP(ss,2),HSP(ss+1,2),nk+1)' linspace(HSP(ss,3),HSP(ss+1,3),nk+1)'];  
  kpoints = [kpoints; p(1:nk,:)];
  list    = [list nk*ss+1];
 
 else
     
  SP = HSP(ss,:);
  p = [linspace(HSP(ss,1),HSP(ss+1,1),nk)' linspace(HSP(ss,2),HSP(ss+1,2),nk)' linspace(HSP(ss,3),HSP(ss+1,3),nk)'];  
  %p = [linspace(SP(1),SP(4),nk)' linspace(SP(2),SP(5),nk)' linspace(SP(3),SP(6),nk)'];  
  kpoints = [kpoints; p(1:nk,:)];
  list    = [list nk*ss];
 
 end
end

kpoints = kpoints*2*pi;
 
nks    = length(kpoints);

end