%%% --- ftn58 of High Order TI with Chiral Hinge --- %%%
%%% --- Science Advances 01 Jun 2018 Vol.4 no.6  --- %%%
%%% ---         DOI: 10.1126/sciadv.aat0346      --- %%%
%%% ---         Hung, Yi-Chun  July/23/2021      --- %%%
%%% ------------------------------------------------ %%%



%% --- Parameters --- %%%
gamma = -0.7;
m1    = 0.6/sqrt(2);

%% --- ftn58 --- %%%
%% m1
ftn58_m1 = [[1 2 (-1i)*m1 0 0 0];[3 4 (-1i)*m1 0 0 0]];

%% gammas
%　gamma_1 
ftn58_gamma_1 = [[1 4 (1i)*(1/2i) 0 1 0];[2 3 (1i)*(1/2i) 0 1 0];...
                 [1 4 (-1)*(1i)*(1/2i) 0 -1 0];[2 3 (-1)*(1i)*(1/2i) 0 -1 0]];
ftn58_gamma_3 = [[1 3 (1i)*(1/2i) 1 0 0];[2 4 (-1i)*(1/2i) 1 0 0];...
                 [1 3 (-1)*(1i)*(1/2i) -1 0 0];[2 4 (1i)*(1/2i) -1 0 0]];
ftn58_gamma_4 = [[1 3 gamma 0 0 0];[1 3 (1/4) 0 0 1];[1 3 (1/2) 1 0 0];...
                                   [1 3 (1/4) 0 0 -1];[1 3 (1/2) -1 0 0]];
ftn58_gamma_4 = [ftn58_gamma_4;[[ftn58_gamma_4(:,1:2)+1],ftn58_gamma_4(:,3:6)]];
ftn58_gamma_2 = [[1 4 gamma 0 0 0];[1 4 (1/4) 0 0 1];[1 4 (1/2) 0 1 0];...
                                   [1 4 (1/4) 0 0 -1];[1 4 (1/2) 0 -1 0]];
ftn58_gamma_2 = [ftn58_gamma_2;...
                [ftn58_gamma_2(:,1)+1,ftn58_gamma_2(:,2)-1,(-1)*ftn58_gamma_2(:,3),ftn58_gamma_2(:,4:6)]];
ftn58_gamma   = [ftn58_gamma_1;ftn58_gamma_2;ftn58_gamma_3;ftn58_gamma_4];

%% ------------------------------------------- %%%
ftn58 = [ftn58_m1;ftn58_gamma];

%%
ftn58_ijdd	 	  = [ftn58(:,1:2) ftn58(:,4:6)];
ftn58_unique_ijdd = unique(ftn58_ijdd,'rows');
[~,Locb] 		  = ismember(ftn58_ijdd,ftn58_unique_ijdd,'rows');
ftn58_unique 	  = [ftn58_unique_ijdd(:,1:2) zeros(size(ftn58_unique_ijdd,1),1) ftn58_unique_ijdd(:,3:5)];
for i=1:length(Locb)
	ftn58_unique(Locb(i),3) = ftn58_unique(Locb(i),3) + ftn58(i,3); 
end
ftn58             = ftn58_unique;
%%

nbond = length(ftn58);
ftn58 = [(1:nbond)' ftn58];
ftn58 = [[4 nbond 0 0 0 0 0];ftn58];

save ftn58.mat ftn58