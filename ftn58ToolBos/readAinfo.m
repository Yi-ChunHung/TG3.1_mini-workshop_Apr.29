load ftn58sparse.mat
Ainfo = struct2table(ftn58sparse.Ainfo);
for j=1:13
	[Ainfo.OrbitIndex{j}]
end
