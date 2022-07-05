% check the maximum hopping distance for bulk TB Hamiltonian
load ftn58sparse.mat

dd = ftn58sparse.dd;
maxHop_a = max(dd(:,1))
maxHop_b = max(dd(:,2))
maxHop_c = max(dd(:,3))