%%% Cut the super cell TBHmftn into Nano-Ribbon %%%
%%% ------------------------------------------- %%%
%%% ------ Hung,YiChun --- Mar.10 (2021) ------ %%%
%%% ------------------------------------------- %%%

clear all

load super_ftn58sparse.mat
load slab_info.mat nv1 nv2 nv3 ribbon_direction

%% --- Initialization 1 --- %%
ij = SPftn58sparse.ij;
tt = SPftn58sparse.tt;
dd = SPftn58sparse.dd;

%%% --- Initialization 2 --- %%
%wcal = ReadInput('input.txt');
%v1   = wcal.vec1(2:end-1);
%sid  = find(isspace(v1));
%vec1 = [str2double(v1(1:sid(1)-1)) str2double(v1(sid(1)+1:sid(2)-1)) ...
%        str2double(v1(sid(2)+1:end))];
%v2   = wcal.vec2(2:end-1);
%sid  = find(isspace(v2));
%vec2 = [str2double(v2(1:sid(1)-1)) str2double(v2(sid(1)+1:sid(2)-1)) ...
%        str2double(v2(sid(2)+1:end))];
%v3   = wcal.vec3(2:end-1);
%sid  = find(isspace(v3));
%vec3 = [str2double(v3(1:sid(1)-1)) str2double(v3(sid(1)+1:sid(2)-1)) ...
%        str2double(v3(sid(2)+1:end))];
%%% ------------------------------------ %%
%%% ---------- dummy evasion ----------- %%
%% is not k (direction that is not periodic) 
%isnk = zeros(1,3);
%if norm(vec1) ~= 1 
%	isnk(1) = 1;
%end
%if norm(vec2) ~= 1
%	isnk(2) = 1;
%end
%if norm(vec3) ~= 1
%	isnk(3) = 1;
%end

for i = 1:3
	if i ~= ribbon_direction 
		indices       = find(dd(:,i)~=0);
		ij(indices,:) = [];
		tt(indices)   = [];
		dd(indices,:) = [];
	end
end

ribdir 					 = zeros(1,3);
ribdir(ribbon_direction) = 1;

NRftn58sparse        = SPftn58sparse;
NRftn58sparse.ij     = ij;
NRftn58sparse.tt     = tt;
NRftn58sparse.dd     = dd(:,ribbon_direction);
%NRftn58sparse.BR1D   = [NRftn58sparse.BR(1,:)*(1/nv1);NRftn58sparse.BR(2,:)*(1/nv2);NRftn58sparse.BR(3,:)*(1/nv3)];
NRftn58sparse.BR1D   = NRftn58sparse.BR(1:3,:).*[(1/nv1) (1/nv2) (1/nv3)];
NRftn58sparse.ribdir = ribdir;
NRftn58sparse        = rmfield(NRftn58sparse,'BR');
NRftn58sparse        = rmfield(NRftn58sparse,'latss');
NRftn58sparse.Orbitps(:,4:6) = NRftn58sparse.Orbitps(:,4:6).*[nv1 nv2 nv3];

for i = 1:size(NRftn58sparse.Ainfo,2)
	position 			  = NRftn58sparse.Ainfo(i).Position;
	position 			  = position.*[nv1 nv2 nv3];
	NRftn58sparse.Ainfo(i).Position   = position;
end

save('Nano_Ribbon_ftn58sparse.mat','NRftn58sparse')
