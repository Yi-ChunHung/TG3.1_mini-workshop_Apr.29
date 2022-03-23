function new_ftn58sparse = ftn58sparse_transform_representation(ftn58sparse,tmode)
% This codes transform the ftn58sparse from Wannier representation 
% into Bloch representation and viece versa
% tmode = 0: Wannier to Blcoh
% tmode = 1: Bloch to Wannier


	%% initialization
	new_ftn58sparse = ftn58sparse;
	ii = ftn58sparse.ij(:,1);
	jj = ftn58sparse.ij(:,2);
	dd = ftn58sparse.dd;
	dd = round(dd,10);
	norb = ftn58sparse.norb;
	isSO = ftn58sparse.isSO;

	%% real procedure
	switch tmode
		case 0
			%if length(find(mod(dd,1))) == 0
			%	fprintf('Please check your ftn58sparse is indeed in Wannier representation \n')
			%	return
			%end
			for atom = 1:length(ftn58sparse.Ainfo)
				orbitindex = ftn58sparse.Ainfo(atom).OrbitIndex;
				if isSO
					orbitindex = [orbitindex,orbitindex+norb/2];
				end
				position   = ftn58sparse.Ainfo(atom).Position;
				v          = find(ismember(ii,orbitindex));
		        dd(v,:)    = dd(v,:) + position;
				w          = find(ismember(jj,orbitindex));
				dd(w,:)    = dd(w,:) - position;
			end
			new_ftn58sparse.dd   = round(dd,10); % to eliminate decimals such as 0 + 1.0e-16 ... etc.
			new_ftn58sparse.type = 'Bloch';
		case 1
			if length(find(mod(dd,1))) ~= 0
				fprintf('Please check your ftn58sparse is indeed in Bloch representation \n')
				return
			end
			for atom = 1:length(ftn58sparse.Ainfo)
				orbitindex = ftn58sparse.Ainfo(atom).OrbitIndex;
				position   = ftn58sparse.Ainfo(atom).Position;
				if isSO
					orbitindex = [orbitindex,orbitindex+norb/2];
				end
				v          = find(ismember(ii,orbitindex));
		        dd(v,:)    = dd(v,:) - position;
				w          = find(ismember(jj,orbitindex));
				dd(w,:)    = dd(w,:) + position;
			end
			new_ftn58sparse.dd   = dd;
			new_ftn58sparse.type = 'Wannier';
		otherwise
			fprintf('tmode is a Boolean !!! \n')
			return
	end
end
