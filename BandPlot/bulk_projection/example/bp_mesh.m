function bp_mesh(ftn58sparse)

%=======================================================================================================
%% --- Read The K-path & Generate the Plotting Data ----- %%%
% --------------------------------------------------------- %
b = strcat('kpt_bp','.labelinfo.dat');
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

BR   = ftn58sparse.BR;
abc  = ftn58sparse.abc;

%% lines to present the length of k-path in BZ
T         = [BR(1,:)*abc(1); BR(2,:)*abc(2); BR(3,:)*abc(3)];
dk_vec    = diff((2*pi*(T\eye(3))*kpt')');
dk_length = zeros(1,length(dk_vec));
for j = 1:length(dk_vec)
    dk_length(j) = norm(dk_vec(j,:));
end
p         = [0,cumsum(dk_length)];
sympt     = p( A.mat(:,1)+1 - cumsum([0;ones(length(A.mat(:,1))-1,1)]) );

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
save('bp_plot.mat','p','kpt','symlb','sympt')
