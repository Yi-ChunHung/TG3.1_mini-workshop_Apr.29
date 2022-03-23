%% --- reread the ij, tt, dd --- %%%
ftn58_temp = [];

j=0;
fid=fopen('slab_ftn58.dat');
while ~feof(fid)
    j=j+1;
    tline=fgetl(fid);

        i=1;
        while tline(i)==' ' %skip blank
            i=i+1;
        end

    temp       = split(tline);
    ij_temp    = str2double(temp(1:2,1))';
    tt_temp    = str2double(temp(3,1)) + 1i*str2double(temp(4,1));
    dd_temp    = str2double(temp(5:7,1))';
    ftn58_temp(j,:) = [ij_temp tt_temp dd_temp];
end
fclose(fid);

ijs = ftn58_temp(:,1:2);
tts = ftn58_temp(:,3);
dds = ftn58_temp(:,4:6);

%% --- Rotating the new unit vectors --- %%
if (isRT)
    BR            = BR*RT'; 
end

%% --- Save as Sftn58sparse format --- %%%
for i=1:ftn58sparse.Nat
    aname = getfield(ftn58sparse.Ainfo, {i}, 'Atom');
end

for i=1:size(superatomps,1)
    atominfo(i).Atom(1,:) = getfield(ftn58sparse.Ainfo, {superatomps(i,5)}, 'Atom');
    atominfo(i).Position(1,1:3) = superatomps(i,9:11);
    atominfo(i).Norb = superatomps(i,4);
    atominfo(i).OrbitIndex = getfield(ftn58sparse.Ainfo, {superatomps(i,5)}, 'OrbitIndex');
    atominfo(i).Orbit(1,:) = getfield(ftn58sparse.Ainfo, {superatomps(i,5)}, 'Orbit');
    atominfo(i).OrbitID    = getfield(ftn58sparse.Ainfo, {superatomps(i,5)}, 'OrbitID');
end

% ob1 = [1:6 7:12 19:24 25:30];
% ob2 = [7:12 13:18 25:30 31:36];
% ttR = 2.5;
% d00 = find(dds(:,1)==0&dds(:,2)==0);
% ii  = ijs(:,1);
% jj  = ijs(:,2);
% for i=1:length(ob1)/3
%     id_ob1{i}  = ob1((i-1)*3+1:i*3);
%     id_ob2{i}  = ob2((i-1)*3+1:i*3);
% end
% for i=1:length(id_ob1)
%     for ib1=1:length(id_ob1{1,1})
%         for ib2=1:length(id_ob2{1,1})
%             hop1 = find(ii(d00)==id_ob1{1,i}(ib1)&jj(d00)==id_ob2{1,i}(ib2));
%             hop2 = find(ii(d00)==id_ob1{1,i}(ib2)&jj(d00)==id_ob2{1,i}(ib1)); 
%             tts(d00(hop1)) = tts(d00(hop1))*ttR;
%             tts(d00(hop2)) = tts(d00(hop2))*ttR;
%         end
%     end
% end

Sftn58sparse.System  = 'Slab_ftn58';
Sftn58sparse.Ainfo   = atominfo;
Sftn58sparse.hkl     = hkl;
Sftn58sparse.abc     = ftn58sparse.abc;
Sftn58sparse.BR3D    = ftn58sparse.BR;
Sftn58sparse.BR2D    = BR;
Sftn58sparse.Nat     = length(atominfo); 
Sftn58sparse.ver     = 'type1';
Sftn58sparse.isSO    = ftn58sparse.isSO;
Sftn58sparse.norb    = size(orbitps,1);
Sftn58sparse.Orbitps = [orbitps(1:end,1:3) orbitps(1:end,4:6) orbitps(1:end,10)];
Sftn58sparse.ij      = ijs;
Sftn58sparse.tt      = tts;
Sftn58sparse.dd      = dds(:,1:2);

save Sftn58sparse.mat Sftn58sparse
