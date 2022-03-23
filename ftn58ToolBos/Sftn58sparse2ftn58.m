function Sftn58=Sftn58sparse2ftn58(Sftn58sparse)

            temp = Sftn58sparse;
            norb = temp.norb;
%            ftn58=[[norb nbond 0 0 0 0 0]; [(1:nbond)' data(:,4:5) data(:,6)+ i*data(:,7) data(:,1:3)]];
            log = find(temp.ij(:,1) > temp.ij(:,2));
	    nbond = size(temp.tt,1);
            Sftn58 = [[norb nbond 0 0 0 0 ]; [(1:nbond)' temp.ij temp.tt temp.dd]];
            Sftn58(log+1,:) = [];
            nbond = nbond-length(log);
            Sftn58(1,2) = nbond;
