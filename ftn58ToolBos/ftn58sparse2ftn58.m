function ftn58=ftn58sparse2ftn58(ftn58sparse)
% <!> the program automatically reduce the ftn58 into upper-triangular form

            temp = ftn58sparse;
            norb=temp.norb;
%            ftn58=[[norb nbond 0 0 0 0 0]; [(1:nbond)' data(:,4:5) data(:,6)+ i*data(:,7) data(:,1:3)]];
            log = find(temp.ij(:,1) > temp.ij(:,2));
	    nbond = size(temp.tt,1);
            ftn58=[[norb nbond 0 0 0 0 0]; [(1:nbond)' temp.ij temp.tt temp.dd]];
            ftn58(log+1,:)=[];
            nbond=nbond-length(log);
            ftn58(1,2)=nbond;
