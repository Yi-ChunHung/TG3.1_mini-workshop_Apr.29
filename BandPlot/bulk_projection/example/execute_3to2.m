
load ftn58sparse.mat 
bp_mesh(ftn58sparse);

BPdirection = 3;    % kx = 1 ; ky = 2 ; kz = 3

amount_of_layer = 100;

PCD_orb_1 = [];
PCD_orb_2 = [];
PCD_orb_3 = [];
PCD_orb_4 = [];

bandplot_concise_ftn58sparse_ex_bp(ftn58sparse, BPdirection, amount_of_layer, PCD_orb_1,PCD_orb_2,PCD_orb_3,PCD_orb_4)

Ef = 0;
Name = 'GaAs';

PlotDispersion(Ef,Name)
