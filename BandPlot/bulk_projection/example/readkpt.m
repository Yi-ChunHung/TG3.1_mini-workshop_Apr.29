function [A]=readkpt(s)
% read from case_band.labelinfo.dat
% to build the k-path
% ----------------------------------------------------------------------------
% label		jth_partition_point			the_k_points
% Gamma     1                           0    0   0
% X         500                         0.5	 0   0
% ... etc.
% ----------------------------------------------------------------------------     
A.mat = [];
position  = 0; % position on the plot axis

j=0;
fid=fopen(s);
while ~feof(fid)
    j=j+1;
    tline=fgetl(fid);

        i=1;
        while tline(i)==' ' %skip blank
            i=i+1;
        end

    temp       = split(tline);
	label      = temp(1,1);
    partition  = cell2mat(temp(2,1));
    partition  = str2double(partition);
    mesh	   = str2double(temp(3:5,1));
    mesh = mesh';
	if j >= 2
		position  = position + norm( abs( mesh - A.mat(j-1,2:4) ));
	end
	A.label(j) = label;
    A.mat(j,:) = [partition,mesh,position];
end
fclose(fid);