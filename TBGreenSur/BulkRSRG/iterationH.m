function [Hs,Hb]=renormH(E,H,T,Ni)
	%iteration process (PRB 31, 5166 eq. 10)
	Hs=H+T/(E-H)*T';
	Hb=H+T/(E-H)*T'+T'/(E-H)*T;
	Tt=T/(E-H)*T;
	Ttp=T'/(E-H)*T';
	for i=1:Ni-1
		Hst=Hs+Tt/(E-Hb)*Ttp;
		Hbt=Hb+Tt/(E-Hb)*Ttp+Ttp/(E-Hb)*Tt;
		Ttt=Tt/(E-Hb)*Tt;
		Ttpt=Ttp/(E-Hb)*Ttp;
		Hs=Hst;
		Hb=Hbt;
		Tt=Ttt;
		Ttp=Ttpt;
	end
end
