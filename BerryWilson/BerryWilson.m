%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates Z2 by Wilson loop method
% base on PRB 84, 075119 (2011)
% 0.01 ver.  author: Tay-Rong Chang  2015/Jan/31
%--------------------------------------------------------------------
% Editor for better ftn58sparse interface: Hung, Yi-Chun 2021/July/27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function BerryWilson

    %% --- Input Arguments --- %%
    wcal = ReadInput('input.txt');
    nVB  = wcal.nVB;
    ib1  = wcal.ib1;
    ib2  = wcal.ib2;
    ib3  = wcal.ib3;
    np   = wcal.ncpu;
    mode = wcal.mode;
    
    if mode ~= 0 && mode ~= 1
        fprintf('Z2: mode=0 ; Chern: mode=1 \n')
        error('There is only two modes !!! \n')
    end
    
    %% --- Load TB-model --- %%
    instruct = load(wcal.ref);
    isSlab   = 0;
    % Sftn58sparse = instruct.Sftn58sparse;
    if isfield(instruct,'ftn58sparse')
        ftn58sparse = instruct.ftn58sparse;
    elseif isfield(instruct,'SPftn58sparse')
        ftn58sparse = instruct.SPftn58sparse;
    elseif isfield(instruct,'Sftn58sparse')
        ftn58sparse = instruct.Sftn58sparse;
        isSlab      = 1;
    end
    
    %% --- Initialization --- %%
    ii   = ftn58sparse.ij(:,1);
    jj   = ftn58sparse.ij(:,2);
    tt   = ftn58sparse.tt;
    dd   = ftn58sparse.dd;
    norb = ftn58sparse.norb;
    theta_Dm = zeros(nVB,ib2,ib3);  
    fff      = zeros(ib2,ib3);
    
    %% --- Generate the Mesh --- %%
    b1=linspace(-0.5,0.5,ib1+1);  % -pi to pi
    b1(end)=[];
    if mode     
        b2=linspace(-0.5,0.5,ib2);  % -pi to pi (BZ for Chern)
        %b2(end)=[]; 
    else
        b2=linspace(0,0.5,ib2);  % 0 to pi (EBZ for Z2)
        %b2(end)=[];
    end
    b3=linspace(-0.5,0.5,ib3);  % 0 to pi
    
    %% --- Opening The Pool --- %%%
    delete(gcp)
    c = parcluster('local');
    c.NumWorkers = np;
    parpool(c, c.NumWorkers);
    
    %% --- Actual Procedure --- %%%
    tic 
    for ii3=1:ib3
       parfor ii2=1:ib2
            %fprintf('tot b2=%d , count=%d \n',ib2,ii2)
    
            F=eye(nVB,nVB);
            f=1;
            for ii1=1:ib1
                %fprintf('Constructing matrix \n')
                %fprintf('tot b1=%d , count=%d \n',ib1,ii1)
                if isSlab 
                    kvec = [b1(ii1),b2(ii2)]*2*pi;
                else   
                    kvec = [b1(ii1),b2(ii2),b3(ii3)]*2*pi;
                end
                H       = full(sparse(ii,jj,exp(1i*dd*kvec').*tt,norb,norb));
                H       = (H+H')/2;
                [EV,E]  = eig(H);
                [~,occ] = sort(diag(E));
                EV_VB   = EV(:,occ(1:nVB));
                %EV_VB  = EV(:,1:nVB);                       % store eigenvector of valence band
                if ii1==1
                    EV1_VB    = EV_VB;
                    EV_pre_VB = EV1_VB;
                else
                    if ii1 == ib1;                           % PRB84,075119 eq. 11: final b1=1st b1
                        EV_VB = EV1_VB;
                    end
                    f = f*det(EV_pre_VB'*EV_VB);             % det(exp(A)) = exp(tr(A))
                    F = F*EV_pre_VB'*EV_VB;                  % PRB84,075119 eq. 10 and 11
                    EV_pre_VB = EV_VB;
                end
            end
            eigval_D = eig(F);                               % PRB84,075119 eq. 11: eigenvalue of D
            theta_Dm(:,ii2,ii3) = sort(angle(eigval_D));     % PRB84,075119 eq. 12      
            fff(ii2,ii3) = angle(f);                                                                  
        end
    end
    toc
    
    save Z2.mat    mode b3 b2 theta_Dm % theta_Dm for all kz slice
    save chern.mat mode b3 b2 fff      % trace of tr(A) on kx,ky plane for all kz slice 
    
    poolobj = gcp('nocreate');
    delete(poolobj);
    
