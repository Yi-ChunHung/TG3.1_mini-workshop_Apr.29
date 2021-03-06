%%%        Program for Calculating Surface Green's Function          %%%
%%% ---------------------------------------------------------------- %%%
%%% Necessary Input files:                                           %%%
%%% Sftn58sparse.mat (includes sufficent # of layers)                %%%
%%% ---------------------------------------------------------------- %%%
%%% Required Functions:                                              %%%
%%% 1) renormH.m                                                     %%%
%%% ---------------------------------------------------------------- %%%
%%% Input Arguments description:                                     %%%
%%% Ni    ==> number of renormalization iteration                    %%%
%%% nn    ==> broadening of calculating Green's Functions            %%%
%%% isPot ==> if true, adding an additional surface layer            %%%
%%% ---------------------------------------------------------------- %%%
%%% Modified by Hans at 25/Jun/2017:                                 %%%
%%% Change ths method for constructing H & T parts by changing the   %%%
%%% orbit order.                                                     %%%
%%% ---------------------------------------------------------------- %%%
function TBSurGreen

    %% --- Input Arguments -- %%
    wcal  = ReadInput('input.txt');
    Ni    = wcal.Ni;
    nn    = wcal.nn;
    Nk    = wcal.Nk; 
    % kp1   = [wcal.kp11 wcal.kp12];
    % kp2   = [wcal.kp21 wcal.kp22];
    % Emu   = [wcal.Emu1,wcal.Emu2];
    Ef    = wcal.Ef;
    nE    = wcal.nE;
    isrev = wcal.isrev;
    np    = wcal.ncpu;
    isPot = wcal.isPot;
    
    k1 = wcal.kp1(2:end-1);
    sid = find(isspace(k1));
    kp1 = [str2double(k1(1:sid(1)-1)) str2double(k1(sid(1)+1:end))];
    
    k2 = wcal.kp2(2:end-1);
    sid = find(isspace(k2));
    kp2 = [str2double(k2(1:sid(1)-1)) str2double(k2(sid(1)+1:end))];
    
    k3 = wcal.kp3(2:end-1);
    sid = find(isspace(k3));
    kp3 = [str2double(k3(1:sid(1)-1)) str2double(k3(sid(1)+1:end))];
    
    em  = wcal.Emu(2:end-1);
    sid = find(isspace(em));
    Emu = [str2double(em(1:sid(1)-1)) str2double(em(sid(1)+1:end))];
    
    
    %%%--------------------------------------------------%%%
    
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
    
    %% --- Initial Setting for E-K dispersion --- %%%
    ef = linspace(Emu(1),Emu(2),nE)+Ef; 
    
    p1 = [linspace(kp1(1),kp2(1),Nk)' linspace(kp1(2),kp2(2),Nk)'];
    p2 = [linspace(kp2(1),kp3(1),Nk)' linspace(kp2(2),kp3(2),Nk)'];
    
    kvec2 =[p1(1:Nk,:);p2(1:Nk,:)];
    %kvec2 =[p1(1:Nk,:)];
    nks = length(kvec2);
    EE  = repmat(ef',1,nks);
    KK  = repmat([1:nks],nE,1);
    %%% ---------------------------------- %%%
    
    %% --- Read Input Information --- %%%
    isSO    = Sftn58sparse.isSO;
    norb    = Sftn58sparse.norb;
    orbitps = Sftn58sparse.Orbitps;
    g_norb  = norb/2;
    ij      = Sftn58sparse.ij;
    dd      = Sftn58sparse.dd(:,1:2);
    tt      = Sftn58sparse.tt;

    %% if isrev => turn the slab upside down
    if isrev
        %% sorting the orbitps in another order
        Orbitps         = Sftn58sparse.Orbitps;
        fh_orbitps      = sortrows(Orbitps(1:(end/2),:),[6],'descend');
        lh_orbitps      = sortrows(Orbitps((end/2 + 1):end,:),[6],'descend');
        reorder_orbitps = [fh_orbitps;lh_orbitps];

        %% map the orbitals according to the new orbitps
        table     = [(1:norb)',reorder_orbitps(:,1)];
        [~,ii_id] = ismember(ij(:,1),table(:,1));
        [~,jj_id] = ismember(ij(:,2),table(:,1));
        ij_ro     = [table(ii_id,2),table(jj_id,2)];
    else
        ij_ro     = ij;
    end
    nij = ij_ro;

    if (isSO)
        Lnorb = g_norb/2;
        upL2  = [Lnorb+1:2*Lnorb];
        dnL1  = [2*Lnorb+1:3*Lnorb];
        for iorb = 1:length(upL2)
            %--- 1st orb
            ib1 = ij_ro(:,1)==upL2(iorb);
            ib2 = ij_ro(:,2)==upL2(iorb);
            nij(ib1,1) = dnL1(iorb);
            nij(ib2,2) = dnL1(iorb);
            %--- 2nd orb
            ib1 = ij_ro(:,1)==dnL1(iorb);
            ib2 = ij_ro(:,2)==dnL1(iorb);
            nij(ib1,1) = upL2(iorb);
            nij(ib2,2) = upL2(iorb);
        end
    end
    
    bulk_id = find(nij(:,1)<(norb/2+1) & nij(:,2)<(norb/2+1));
    hopp_id = find(nij(:,1)<(norb/2+1) & nij(:,2)>(norb/2));
    ii_b    = nij(bulk_id,1);
    jj_b    = nij(bulk_id,2);
    ii_s    = nij(hopp_id,1);
    jj_s    = nij(hopp_id,2)-(norb/2);
    
    dd_b    = dd(bulk_id,:);
    tt_b    = tt(bulk_id,:);
    dd_s    = dd(hopp_id,:);
    tt_s    = tt(hopp_id,:);
    %%% --------------------------------------------------- %%%
    
    %% --- Create a New surface Hamiltonian (with different potential) --- %
    if isPot == 1
        topsur = find(orbitps(:,6)>wcal.dis&orbitps(:,1)<=(norb/2));
        topsur = [topsur;topsur+(max(ii_b)/2)];
        nH_id  = find(dd_b(:,1)==0&dd_b(:,2)==0&ii_b(:)==jj_b(:));
        nH_top_id = [];
        for ii=1:length(topsur)
            tmp_top_id = find(ii_b(nH_id(:))==topsur(ii));
            nH_top_id  = [nH_top_id;tmp_top_id];
        end
        nH_tt = tt_b;
        nH_tt(nH_id(nH_top_id),:) = nH_tt(nH_id(nH_top_id),:) + wcal.PeV;
    end
    %%% ------------------------------------------------------- %%%
    
    %% --- Pauli Matrix %%
    sx0 = [0 1;1 0];
    sy0 = [0 -1i;1i 0];
    sz0 = [1 0;0 -1];
    Sx = kron(sx0,eye(g_norb/2));
    Sy = kron(sy0,eye(g_norb/2));
    Sz = kron(sz0,eye(g_norb/2));
    
    %% --- Actual Procedure --- %%%
    c = parcluster('local');
    c.NumWorkers = np;
    parpool(c, c.NumWorkers);
    
    Ab  = zeros(length(ef),nks);
    As  = zeros(length(ef),nks);
    SSx = zeros(length(ef),nks);
    SSy = zeros(length(ef),nks);
    SSz = zeros(length(ef),nks);
    A00 = zeros(length(ef),nks);
    SSx0 = zeros(length(ef),nks);
    SSy0 = zeros(length(ef),nks);
    SSz0 = zeros(length(ef),nks);
    II  = eye(g_norb,g_norb);
    
    tic
    for mm=1:length(ef)
        if (mod(mm,10)==0)
            fprintf('# of E points = %4i/%i\n',mm,length(ef));
        end
        runw = ef(1,mm);
        Abt  = zeros(length(kvec2),1);
        Ast  = zeros(length(kvec2),1);
        Sxt  = zeros(length(kvec2),1);
        Syt  = zeros(length(kvec2),1);
        Szt  = zeros(length(kvec2),1);
        
        A00t = [];
        Sxt0 = zeros(length(kvec2),1);
        Syt0 = zeros(length(kvec2),1);
        Szt0 = zeros(length(kvec2),1);
        parfor ik=1:length(kvec2)
     %       fprintf('# of k points = %4i/%i\n',ik,length(kvec2));
            kcolumnvec = kvec2(ik,:)'*2*pi;
            E  = (runw+1i*nn)*II;
            H  = full(sparse(ii_b,jj_b,exp(1i*dd_b*kcolumnvec).*tt_b,g_norb,g_norb));
            H  = (H+H')/2;
            T  = full(sparse(ii_s,jj_s,exp(1i*dd_s*kcolumnvec).*tt_s,g_norb,g_norb));
               
            %%% Green't function (renomalization procedure)
            [Hs,Hb] = renormH(E,H,T,Ni);
            Gb      = (E-Hb)\II;
            Gs      = (E-Hs)\II;
            Abt(ik) = -imag(trace(Gb));
            Ast(ik) = -imag(trace(Gs));
            
            Sxt(ik) = -imag(trace(Sx*Gs));
            Syt(ik) = -imag(trace(Sy*Gs));
            Szt(ik) = -imag(trace(Sz*Gs));
            
            if isPot==1
                nH       = full(sparse(ii_b,jj_b,exp(1i*dd_b*kcolumnvec).*nH_tt,g_norb,g_norb)); 
                nH       = (nH+nH')/2;
                nHS      = T*Gs*T';
                GnH      = (E-nH-nHS)\II;
                A00t(ik) = -imag(trace(GnH));
                Sxt0(ik) = -imag(trace(Sx*GnH));
                Syt0(ik) = -imag(trace(Sy*GnH));
                Szt0(ik) = -imag(trace(Sz*GnH));
            end
        end
        
        As(mm,:) = Ast;
        Ab(mm,:) = Abt;
        SSx(mm,:) = Sxt;
        SSy(mm,:) = Syt;
        SSz(mm,:) = Szt;
        
        if isPot==1
            A00(mm,:) = A00t;
            SSx0(mm,:) = Sxt0;
            SSy0(mm,:) = Syt0;
            SSz0(mm,:) = Szt0;
        end
    end
    save('Ek.mat','As','Ab','A00','SSx','SSy','SSz','SSx0','SSy0','SSz0','KK','EE','ef','Ef','Nk','-v7.3') 
    toc
     
    poolobj = gcp('nocreate');
    delete(poolobj);
    
    %% Plotting %%
    % ncolor = 1e4;
    % map    = [linspace(1,0,ncolor)' linspace(1,0,ncolor)' linspace(1,1,ncolor)'];
    % 
    % 
    % if isPot==1
    %    figure,pcolor(KK,EE,A00(:,:)/pi),shading interp;
    % end
    % 
    % if Ktype==1
    %     figure,pcolor(KK,EE-Ef,As(:,:)/pi),shading interp;
    %     % axis equal;
    %     colormap(map);
    %     caxis([0 10]);
    %     save Ek_010_s1.0.mat As Ab A00 KK EE ef
    % else
    %     figure,pcolor(X,Y,As(:,:)/pi),shading interp;
    %     % axis equal;
    %     colormap(map);
    %     caxis([0 200]);
    %     save ES.mat As Ab A00 X Y ef
    % end