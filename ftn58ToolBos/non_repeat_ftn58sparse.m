function unique_ftn58sparse = non_repeat_ftn58sparse(ftn58sparse)
	% find the ftn58sparse with non-repeated [ij,dd]
	% 'ispara' determines whether using parallel computation or not
	
	ftn58_ijdd	 	  = [ftn58sparse.ij ftn58sparse.dd];
	ftn58_unique_ijdd = unique(ftn58_ijdd,'rows');
	[~,Locb] 		  = ismember(ftn58_ijdd,ftn58_unique_ijdd,'rows');

	%if ispara
	%	parfor i=1:size(ftn58_unique_ijdd,1)
	%		Lia               = ismember(Locb,i); 
	%		ftn58_unique(i,6) = sum(ftn58sparse.tt(find(Lia))); 
	%	end
	%else
	%	for i=1:length(Locb)
	%		ftn58_unique(Locb(i),6) = ftn58_unique(Locb(i),6) + ftn58sparse.tt(i); 
	%	end
	%end
	
	unique_ftn58sparse    = ftn58sparse;
	unique_ftn58sparse.ij = ftn58_unique_ijdd(:,1:2);
	unique_ftn58sparse.tt = accumarray(Locb,ftn58sparse.tt)
	unique_ftn58sparse.dd = ftn58_unique_ijdd(:,3:5);
	
	
	
	