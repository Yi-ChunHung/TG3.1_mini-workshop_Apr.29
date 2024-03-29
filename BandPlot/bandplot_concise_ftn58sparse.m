function bandplot_concise_ftn58sparse(ftn58sparse, kpt_filename, PCD_orb_1,PCD_orb_2,PCD_orb_3,PCD_orb_4)
% kpt_filename : the name of the file '.labelinfo.dat'
% more detail information please see the help of 'readkpt.m'
% PCD_orb : orbitals that you want to calculate their partial charge distributions
% e.g. PCD_orb = [1 2 3];
% PCD_orb = [] to NOT have the calculation about partial charge distributions
% add PCD_orb_2 to plot another partial charge distribution in case ... etc.


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
%=======================================================================================================
%% --- Read The K-path & Generate the Plotting Data ----- %%%
% --------------------------------------------------------- %
b = strcat(kpt_filename,'.labelinfo.dat');
A = readkpt(b);
%------------------------------------------
% A.label(j)= label;
% A.mat(j,:)= [partition;mesh;position];
%------------------------------------------
n = size(A.mat,1);
symlb = A.label;
for j = 1:n-1
	%p_temp   = linspace(A.mat(j,5),A.mat(j+1,5),A.mat(j+1,1)-A.mat(j,1));                  % lines that do not present the length of k-path in BZ
	kpt_temp = [linspace(A.mat(j,2),A.mat(j+1,2),A.mat(j+1,1)-A.mat(j,1))' ...
				linspace(A.mat(j,3),A.mat(j+1,3),A.mat(j+1,1)-A.mat(j,1))' ...
				linspace(A.mat(j,4),A.mat(j+1,4),A.mat(j+1,1)-A.mat(j,1))' ];				
	if j == 1
	%	p     = p_temp;
		kpt   = kpt_temp;
	%	sympt = A.mat(j,5);
	else
	%	p     = [p(1:end-1),p_temp];
		kpt   = [kpt(1:end-1,:);kpt_temp];
	%	sympt = [sympt, A.mat(j,5)];
	end
end

%% lines to present the length of k-path in BZ
T         = diag(abc)*BR; %[BR(1,:)*abc(1); BR(2,:)*abc(2); BR(3,:)*abc(3)];
dk_vec    = diff((2*pi*(T\eye(3))*kpt')');
dk_length = zeros(1,length(dk_vec));
for j = 1:length(dk_vec)
    dk_length(j) = norm(dk_vec(j,:));
end
p         = [0,cumsum(dk_length)];
sympt     = p( A.mat(:,1)+1 - cumsum([0;ones(length(A.mat(:,1))-1,1)]) );

%disp(p)
%disp(kpt)
% --------------------------------------------------------- %
%% --- Paragraph Output: 'p', 'kpt', 'sympt', 'symlb' --- %%%
% --------------------------------------------------------- %
% p     : the list records the position of each k-points on the axis of the dispersion plot
% kpt   : the list records the k-points along the choosen k-path
% sympt : the list records the position of the high symmetry points on the axis of the dispersion plot
% symlb : the list records the corresponding characters of the high symmetry points
% ------------------------------------------------------------------------------------------------------
% p     = [ p1, p2, ... ];
% kpt   = [ [kpt1.x kpt1.y kpt1.z]; [kpt2.x kpt2.y kpt2.z]; ... ];
% sympt = [ sympt1, sympt2, ... ];
% symlb = [ 'Gamma', 'X', ... ];
%=======================================================================================================
%% --- Calculate the Eigenvalues, Eigenvectors and the Partial Charge Distriution ---- %%%
% -------------------------------------------------------------------------------------- %
% ---------------------- Preallocation ---------------------- % 
% ----------------------------------------------------------- % 
nks       = size(kpt,1);
kpoints   = 2*pi*kpt;
Ek        = zeros(nks,norb);
PCD_1 	  = zeros(nks,norb);
PCD_2     = zeros(nks,norb);
PCD_3     = zeros(nks,norb);
PCD_4     = zeros(nks,norb);

% ----------------------------------------------------------- %
% ---------------------- Calculation ------------------------ % 
% ----------------------------------------------------------- %
tic
fprintf('Start calculating ... \n')
parfor ik=1:nks
    HH_sparse = sparse(ii,jj,exp(1i*dd*kpoints(ik,:)').*tt,norb,norb);
    HH        = full(HH_sparse);
    HH        = (HH+HH')/2;
    [vec, Etemp] = eig(HH);
     
	% Record the calculation 
    % Ham{ik,1}    = HH;
    % eigvec{ik,1} = vec;
    
    Ek(ik,:)  = diag(Etemp);
	% Calculate the partial charge distributions
    if length(PCD_orb_1) ~= 0
		PCD_1(ik,:) = vecnorm(vec(PCD_orb_1,:)).^2;
    end 
    if length(PCD_orb_2) ~= 0
		PCD_2(ik,:) = vecnorm(vec(PCD_orb_2,:)).^2;
    end
    if length(PCD_orb_3) ~= 0
        PCD_3(ik,:) = vecnorm(vec(PCD_orb_3,:)).^2;
    end
    if length(PCD_orb_4) ~= 0
        PCD_4(ik,:) = vecnorm(vec(PCD_orb_4,:)).^2;
    end

end
toc

save Ekplot.mat norb p symlb sympt Ek PCD_1 PCD_2 PCD_3 PCD_4 PCD_orb_1 PCD_orb_2 PCD_orb_3 PCD_orb_4 

%close the pool
delete(gcp);
return
