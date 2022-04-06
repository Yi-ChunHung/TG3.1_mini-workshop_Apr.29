function TBGreenSur_bulkRSRG()

    % ----- Input parameter -----%
    wcal    = ReadInput('input.txt');
    Ni      = wcal.Ni;                    % # of interation (layer); layer=2^Ni
    ef      = wcal.ef;
    Ne      = wcal.Ne;
    np      = wcal.np;
    nn      = wcal.nn;                  % energy broadening
    kc      = wcal.kc;
    surface = wcal.surface; % surface x=1, y=2, z=3
    hop_d   = wcal.hop_d;   % hopping direction. positive=1(coord small surface), negeative=2;
    %e.g hop_d = 1 means that the lower atoms hop to highet atoms , s.t the
    %R>0, and using these T to calculate semi-infinite Green , we have the 
    %semi-infinite Green function which terminated at the lower surface
    %------------- if hop_d=2(semi-infinite at this surface)
    %   BulkBulkBulk
    %   BulkBulkBulk
    %   BulkBulkBulk
    %   BulkBulkBulk
    %------------- if hop_d=1(semi-infinite at this surface)
    
    %% --- Load TB-model --- %%s
    instruct = load(wcal.file);
    % Sftn58sparse = instruct.Sftn58sparse;
    if isfield(instruct,'ftn58sparse')
        Sftn58sparse = instruct.Sftn58sparse;
    elseif isfield(instruct,'SPftn58sparse')
        Sftn58sparse = instruct.SPftn58sparse;
    else
        Sftn58sparse = instruct.Sftn58sparse;
    end
    ftn58 = [(1:length(ftn58sparse.tt))' ftn58sparse.ij ftn58sparse.tt ftn58sparse.dd];
    ftn58 = [[ftn58sparse.norb length(ftn58sparse.tt) 0 0 0 0 0];ftn58];

    % ================================================================== %
    norb           = ftn58(1,1);     % # of orbital
    %--- Choose the PDOS orbitals (OPTIONAL) --- %
    PCD_orb        = 1:norb;
    v              = ones(length(PCD_orb),1);
    PCD_exp        = full(sparse(PCD_orb,PCD_orb,v,norb,norb));
    % ================================================================== %

    k1  = wcal.kp1(2:end-1);
    sid = find(isspace(k1));
    kp1 = [str2double(k1(1:sid(1)-1)) str2double(k1(sid(1)+1:end))];

    k2  = wcal.kp2(2:end-1);
    sid = find(isspace(k2));
    kp2 = [str2double(k2(1:sid(1)-1)) str2double(k2(sid(1)+1:end))];
    
    k3  = wcal.kp3(2:end-1);
    sid = find(isspace(k3));
    kp3 = [str2double(k3(1:sid(1)-1)) str2double(k3(sid(1)+1:end))];

    em  = wcal.Emu(2:end-1);
    sid = find(isspace(em));
    Emu = [str2double(em(1:sid(1)-1)) str2double(em(sid(1)+1:end))];

    tic
    w     = ef+linspace(Emu(1), Emu(2),Ne);    % energy window
    kx    = [linspace(kp1(1),kp2(1),kc),linspace(kp2(1),kp3(1),kc)];
    ky    = [linspace(kp1(2),kp2(2),kc),linspace(kp2(2),kp3(2),kc)];
    kpath = [kx' ky' zeros(2*kc,1)];
    nk    = length(kpath);        % # of k

    if surface==1;
        surf_d = 5;  % ftn58 x
    end
    if surface==2;
        surf_d = 6;  % ftn58 y
    end
    if surface==3;
        surf_d = 7;  % ftn58 z
    end
    
    %% longest hopping in bulk
    long_hop_p = max(max(ftn58(:,surf_d))); % positive hop 
    long_hop_n = min(min(ftn58(:,surf_d))); % negative hop
    
    % brick 
    %brick_ind=find(ftn58(:,surf_d)<long_hop_p & ftn58(:,surf_d)>long_hop_n); % search element of H1 in PRB 31, 5166 eq.12
    
    % the ftn58 of the 
    if hop_d==1;
        brick_ind = find(ftn58(:,surf_d)<long_hop_p & ftn58(:,surf_d)>=0);
    else
        brick_ind = find(ftn58(:,surf_d)>long_hop_n & ftn58(:,surf_d)<=0);
    end
    
    ftn58_brick = ftn58(brick_ind(:),:);
    %n_element_brick=size(brick_ind);         % # of elements
    %ftn58_brick=zeros(n_element_brick(1),7); % memory
    %for iib=1:n_element_brick(1)             % element of brick (ftn58 fromat)
    %    ftn58_brick(iib,:)=ftn58(brick_ind(iib),:);
    %end
    
    % positive hopping (brick to brick)
    if hop_d==1;
        hop_ind_p   = find(ftn58(:,surf_d)>0);    % search element of T1 in PRB 31, 5166 eq.12
        ftn58_hop_p = ftn58(hop_ind_p(:),:);
        %n_element_hop_p=size(hop_ind_p);      % # of elements 
        %ftn58_hop_p=zeros(n_element_hop_p(1),7); % memory 
        %for iih=1:n_element_hop_p(1)             % element of hopping (ftn58 fromat)  
        %    ftn58_hop_p(iih,:)=ftn58(hop_ind_p(iih),:);
        %end
        ftn58_hop = [ftn58(1,:); ftn58_hop_p];
    else
    % negative hopping (brick to brick)
        hop_ind_n   = find(ftn58(:,surf_d)<0);    
        ftn58_hop_n = ftn58(hop_ind_n(:),:);
        %n_element_hop_n=size(hop_ind_n);
        %ftn58_hop_n=zeros(n_element_hop_n(1),7);
        %for iih=1:n_element_hop_n(1)
        %    ftn58_hop_n(iih,:)=ftn58(hop_ind_n(iih),:);
        %end
        ftn58_hop = [ftn58(1,:); ftn58_hop_n];
    end
    
    % =================================================================== %
    % ------- Construct Matrix -------%
    fprintf('Constructing matrix \n');
    for ik=1:length(kpath)
    % parfor ik=1:length(kpath)
        fprintf('matrix total k=%d, now=%d \n',length(kpath),ik);
        kcolumnvec=kpath(ik,:)'*2*pi;
        
        % ----- bulk part -----%
        % ----- transfor ftn58 to matirx (element) ----%
        
        % ---- diagonal element H ----%
        Hk{ik}=full(sparse(ftn58_brick(2:end,2),...
            ftn58_brick(2:end,3),...
            exp(i*ftn58_brick(2:end,5:7)*kcolumnvec).*ftn58_brick(2:end,4),...
            norb,norb));
        Hk{ik}=(Hk{ik}+Hk{ik}')/2;   % keep Hermination
        
        %     % --- off-diagonal term T ---%
        Tk{ik}=full(sparse(ftn58_hop(2:end,2),...
            ftn58_hop(2:end,3),...
            exp(i*ftn58_hop(2:end,5:7)*kcolumnvec).*ftn58_hop(2:end,4),...
            norb,norb));
        
        %     % ---- surface part ----%
        %     % ---- diagonal element H ----%
        Hs0{ik}=full(sparse(ftn58_brick(2:end,2),...
            ftn58_brick(2:end,3),...
            exp(i*ftn58_brick(2:end,5:7)*kcolumnvec).*ftn58_brick(2:end,4),...
            norb,norb));
        Hs0{ik}=(Hs0{ik}+Hs0{ik}')/2;
        % --- off-diagonal term T ---%
        Ts0{ik}=full(sparse(ftn58_hop(2:end,2),...
            ftn58_hop(2:end,3),...
            exp(i*ftn58_hop(2:end,5:7)*kcolumnvec).*ftn58_hop(2:end,4),...
            norb,norb));
        
    end
    
    % --- Memory preallocation ---- %
    II    = eye(norb);
    sx    = kron([0 1;1 0],eye(norb/2));
    sy    = kron([0 -1i;1i 0],eye(norb/2));
    sz    = kron([1 0;0 -1],eye(norb/2));
    Ab    = zeros(length(w),length(kpath));
    As    = zeros(length(w),length(kpath));
    As_sx = zeros(length(w),length(kpath));
    As_sz = zeros(length(w),length(kpath));
    As_sy = zeros(length(w),length(kpath));
    %---------------------%
    
    %% --- Actual Procedure --- %%%
    c = parcluster('local');
    c.NumWorkers = np;
    parpool(c, c.NumWorkers);
    
    %------- RSRG and spectral weight ---- %
    fprintf('Start scanning ...')
    %for iw=1:length(w)    % scan energy
    parfor iw=1:length(w)    % scan energy
        %fprintf('total w=%d, now=%d \n',length(w),iw);
        runw=w(iw);
        Abt=Ab(iw,:);
        Ast=As(iw,:);
        Ast_sx=As_sx(iw,:);
        Ast_sy=As_sy(iw,:);
        Ast_sz=As_sz(iw,:);
        
        for ik=1:length(kpath)    % scan k
            E=(runw+1i*nn)*II;    % w in geen function
            H=full(Hk{ik});
            T=full(Tk{ik});
            H0=full(Hs0{ik});
            T0=full(Ts0{ik});
            [Hs,Hb]=iterationH(E,H,T,Ni);  % PRB 31, 5166 eq.10
            
            %Gb=inv(E-Hb);  % bulk  % PRB 31, 5166 eq.5
            %Gs=inv(E-Hs);  % surface
            Gb=(E-Hb)\II;
            Gs=(E-Hs)\II;
            
            %Gs0=inv(E-H0-T0*Gs*T0');    % cover surface % PRB 31, 5166 eq.13
            Gs0 = (E-H0-T0*Gs*T0')\II;
            
            Abt(ik)=-imag(trace(Gb*PCD_exp));   % spectral weight of bulk
            Ast(ik)=-imag(trace(Gs0*PCD_exp));  % spectral weight of surface (sun over all of total orbital)
                    
            Gs0sx=sx*Gs0;
            Gs0sy=sy*Gs0;
            Gs0sz=sz*Gs0;
            Ast_sx(ik)=-imag(trace(Gs0sx*PCD_exp));
            Ast_sy(ik)=-imag(trace(Gs0sy*PCD_exp));
            Ast_sz(ik)=-imag(trace(Gs0sz*PCD_exp));
        end
        
        Ab(iw,:)=Abt;
        As(iw,:)=Ast;
        As_sx(iw,:)=Ast_sx;
        As_sy(iw,:)=Ast_sy;
        As_sz(iw,:)=Ast_sz;
        
    end
    % ----------------------------------- %
    
    save Ek.mat  w nk ef Ab As As_sx As_sy As_sz
    % save G_matrix.mat -v7.3 Gb
    toc
    
    poolobj = gcp('nocreate');
    delete(poolobj);
return
