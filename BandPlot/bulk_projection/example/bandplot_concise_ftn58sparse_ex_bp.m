function bandplot_concise_ftn58sparse_ex_bp(ftn58sparse, BPdirection , amount_of_layer, PCD_orb_1,PCD_orb_2,PCD_orb_3,PCD_orb_4)
% kpt_filename : the name of the file '.labelinfo.dat'
% more detail information please see the help of 'readkpt.m'
% PCD_orb : orbitals that you want to calculate their partial charge distributions
% e.g. PCD_orb = [1 2 3];
% PCD_orb = [] to NOT have the calculation about partial charge distributions
% add PCD_orb_2 to plot another partial charge distribution in case ... etc.
% amount_of_layer determine the number of partitions along the projected directions 

if length(find(ismember(BPdirection,[1 2 3])))~=1
	fprintf('Check your Bulk Projection direction !!! \n')
	fprintf('There should only be one direction !!! \n')
	return
end

%% Open a pool for parallel computation %%%
% --------------------------------------- %
delete(gcp);
pool=parpool('local');
pool.IdleTimeout=1800; %minutes
%=======================================================================================================
%% ----------------- Initialization --------------------- %%%
% --------------------------------------------------------- %
norb = ftn58sparse.norb;
ii   = ftn58sparse.ij(:,1);
jj   = ftn58sparse.ij(:,2);
dd   = ftn58sparse.dd;
tt   = ftn58sparse.tt;
BR   = ftn58sparse.BR;
abc  = ftn58sparse.abc;
Sz   = [1 0;0 -1];
Sy   = [0 -1i;1i 0];
Sx   = [0 1;1 0];
%=======================================================================================================
load bp_plot.mat
%=======================================================================================================
nks   = size(kpt,1);
Ek    = zeros(nks,norb,amount_of_layer);
PCD_1 = zeros(nks,norb,amount_of_layer);
PCD_2 = zeros(nks,norb,amount_of_layer);
PCD_3 = zeros(nks,norb,amount_of_layer);
PCD_4 = zeros(nks,norb,amount_of_layer);
for layer=1:amount_of_layer
	%=======================================================================================================
	%% --- Calculate the Eigenvalues, Eigenvectors and the Partial Charge Distriution ---- %%%
	% -------------------------------------------------------------------------------------- %
	% ---------------------- Preallocation ---------------------- % 	
	% ----------------------------------------------------------- % 
	if BPdirection == 3
		kpt_temp  = [kpt(:,1), kpt(:,2), (layer-1)*(2*pi/(amount_of_layer+1))*ones(nks,1) ];
	elseif BPdirection == 2
		kpt_temp  = [kpt(:,1), (layer-1)*(2*pi/(amount_of_layer+1))*ones(nks,1), kpt(:,2) ];
	elseif BPdirection == 1
		kpt_temp  = [(layer-1)*(2*pi/(amount_of_layer+1))*ones(nks,1), kpt(:,1), kpt(:,2)];
	end
	kpoints   = 2*pi*kpt_temp;
	PCD_exp_1 = zeros(norb,norb);
	PCD_exp_2 = zeros(nks,norb);
	PCD_exp_3 = zeros(nks,norb);
	PCD_exp_4 = zeros(nks,norb);
	% ----------------------------------------------------------- %
	% --- Expectation Value for the partial charge distribution - % 
	% --- Using block diagonalized Identity matrix -------------- %
	% ----------------------------------------------------------- %
	if length(PCD_orb_1) ~= 0
		v              = ones(length(PCD_orb_1),1);
		PCD_exp_1      = full(sparse(PCD_orb_1,PCD_orb_1,v,norb,norb)); 
	end
	if length(PCD_orb_2) ~= 0
		v              = ones(length(PCD_orb_2),1);
		PCD_exp_2      = full(sparse(PCD_orb_2,PCD_orb_2,v,norb,norb)); 
	end
	if length(PCD_orb_3) ~= 0
		v              = ones(length(PCD_orb_3),1);
		PCD_exp_3      = full(sparse(PCD_orb_3,PCD_orb_3,v,norb,norb));
	end
	if length(PCD_orb_4) ~= 0
		v              = ones(length(PCD_orb_4),1);
		PCD_exp_4      = full(sparse(PCD_orb_4,PCD_orb_4,v,norb,norb));
	end           
	% ----------------------------------------------------------- %
	% ---------------------- Calculation ------------------------ % 
	% ----------------------------------------------------------- %
	%tic
	%fprintf('Start calculating ... \n')
	layer = layer;
	fprintf('%f/%f \n',layer,amount_of_layer)
	parfor ik=1:nks
	    Hsparse   = sparse(ii,jj,exp(1i*dd*kpoints(ik,:)').*tt,norb,norb);
	    HH        = full(Hsparse);
	    HH        = (HH+HH')/2;
	    [vec, Etemp] = eig(HH);
	     
		% Record the calculation 
	    % Ham{ik,1}    = HH;
	    % eigvec{ik,1} = vec;
	    
	    Ek(ik,:,layer)  = diag(Etemp);
		% Calculate the partial charge distributions
	    if length(PCD_orb_1) ~= 0
			PCD_1(ik,:,layer) = diag(vec' * PCD_exp_1 * vec);
	    end 
	    if length(PCD_orb_2) ~= 0
			PCD_2(ik,:,layer) = diag(vec' * PCD_exp_2 * vec);
	    end
	    if length(PCD_orb_3) ~= 0
			PCD_3(ik,:,layer) = diag(vec' * PCD_exp_3 * vec);
	    end
	    if length(PCD_orb_4) ~= 0
			PCD_4(ik,:,layer) = diag(vec' * PCD_exp_4 * vec);
	    end
	end
	Ek = Ek;
	%toc
end

save Ekplot.mat amount_of_layer norb p Ek PCD_1 PCD_2 PCD_3 PCD_4 PCD_orb_1 PCD_orb_2 PCD_orb_3 PCD_orb_4 symlb sympt 

%close the pool
delete(gcp);
	
return