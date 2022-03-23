function [H,trsl] = RS_Hamiltonian(ftn58sparse,np)
%% the function tidy up the matrix elements of a ftn58sparse for each translation

    if np < 1
        fprintf('number of cpus must > 0 !!! \n')
        return
    end

    %% initialization %%
    trsl = unique(new_ftn58sparse.dd,'rows');
    nlen = length(trsl)
    dd   = new_ftn58sparse.dd;
    ii   = new_ftn58sparse.ij(:,1);
    jj   = new_ftn58sparse.ij(:,2);
    tt   = new_ftn58sparse.tt;
    norb = new_ftn58sparse.norb;

    %% real procedure
    H    = cell(1,nlen);
    switch np
    case 1
        for j=1:nlen
            ordnj = find(ismember(dd,trsl(j,:),'rows'));
            H{j}  = full(sparse(ii(ordnj),jj(ordnj),tt(ordnj),norb,norb));
        end
        H = reshape(cell2mat(H),norb,norb,nlen);
    otherwise
        c = parcluster('local');
        c.NumWorkers = np;
        parpool(c, c.NumWorkers);

        parfor j=1:nlen
            ordnj = find(ismember(dd,trsl(j,:),'rows'));
            H{j}  = full(sparse(ii(ordnj),jj(ordnj),tt(ordnj),norb,norb));
        end
        H = reshape(cell2mat(H),norb,norb,nlen);

        poolobj = gcp('nocreate');
        delete(poolobj);
    end
end
