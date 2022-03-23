function check_onsite(ftn58sparse)
% check the onsite of a ftn58saprse

    ij = ftn58sparse.ij;
    tt = ftn58sparse.tt;
    dd = ftn58sparse.dd;
    ftn58 = [ij tt dd];

    on_id  = find(ij(:,1)==ij(:,2) & dd(:,1)==0 & dd(:,2)==0 & dd(:,3)==0);
    onsite = ftn58(on_id,2:3);

    T = array2table(onsite,'VariableNames',{'Spin_Orbitals','Onsite_Energy'})
end