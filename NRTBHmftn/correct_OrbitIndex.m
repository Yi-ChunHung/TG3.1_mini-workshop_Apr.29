function Sftn58sparse = correct_OrbitIndex(Sftn58sparse) 
% The codes correct the Ainfo.OrbitIndex in Sftn58sparse 

    Orbitps = Sftn58sparse.Orbitps;
    isSO    = Sftn58sparse.isSO;
    norb    = Sftn58sparse.norb;
    for atom = 1:Sftn58sparse.Nat
        indices = find(Orbitps(:,3)==atom);
        OrbitIndex = Orbitps(indices,1)';
        if isSO
            OrbitIndex(OrbitIndex>(norb/2))=[];          % neglect the other spin for convention
        end
        Sftn58sparse.Ainfo(atom).OrbitIndex = OrbitIndex;
    end