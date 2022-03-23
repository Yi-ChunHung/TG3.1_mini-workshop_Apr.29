
load ftn58sparse.mat 

kpt_filename = 'kpt';

PCD_orb_1 = [];     % partial charge distribution 
PCD_orb_2 = [];     % input orbit index form Ainfo
PCD_orb_3 = [];     % with a single index or an array of indeces
PCD_orb_4 = [];

bandplot_concise_ftn58sparse(ftn58sparse, kpt_filename, PCD_orb_1,PCD_orb_2,PCD_orb_3,PCD_orb_4)

Ef = 0;
Name = 'GaAs';

PlotDispersion(Ef,Name)
