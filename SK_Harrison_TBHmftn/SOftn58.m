function [newftn58, ftn58SO]=SOftn58(iob,ftn58,SOp,ipp,SOd,ipd)
%%% orbital order:
%%% ipp is the index of px, and for the order px,py,pz
%%% ipd is the index of dz2, and for the order dz2,dxz,dyz,dxy,dx2y2 
%%% ----------------------------------------------------------------- %%%
%%%    iob==0, both p and d orbitals are considered
%%%    iob==1, only p orbitals are considered 
%%%    iob==2, only d orbitals are consideerd
%%% ----------------------------------------------------------------- %%%
norb=ftn58(1,1);
nbond=ftn58(1,2);

ftn58hole=ftn58(2:end,:);
ftn58hole(:,2:3)=ftn58hole(:,2:3)+norb;

newftn58=[ftn58;ftn58hole];
newftn58(1,1)=norb*2;
newftn58(1,2)=nbond*2;
newftn58(2:end,1)=(1:newftn58(1,2))';

nw=size(ftn58,2);

%%%--- For p orbitals only
if nargin <= 5 && iob==1
    
    ngap_p=length(ipp);
    if length(SOp) == 1
        SOp=repmat(SOp,ngap_p,1);
    end
    
    ftn58temp=zeros([ngap_p*6 nw]);
    for ii=1:ngap_p
        ftn58temp((ii*6-5):(ii*6),2:4)= ...
             [ipp(ii) ipp(ii)+1 -1i*SOp(ii)/3            % x(u) y(u)
              ipp(ii) ipp(ii)+norb+2 SOp(ii)/3           % x(u) z(d)
              ipp(ii)+1 ipp(ii)+norb+2 -1i*SOp(ii)/3     % y(u) z(d)
              ipp(ii)+2 ipp(ii)+norb -SOp(ii)/3          % z(u) x(d)
              ipp(ii)+2 ipp(ii)+norb+1 1i*SOp(ii)/3      % z(u) y(d)
              ipp(ii)+norb ipp(ii)+norb+1 1i*SOp(ii)/3]; % x(d) y(d) 
    end

    ftn58SO         = [zeros([1 nw]); ftn58temp];
    ftn58SO(1,1)    = norb*2;
    ftn58SO(1,2)    = ngap_p*6;
    ftn58SO(2:end,1)= (1:ftn58SO(1,2))';
    
%%%--- For d orbitals only
elseif nargin <= 5 && iob==2
    
    ngap_d=length(ipd);
    if length(SOd) == 1
        SOd=repmat(SOd,ngap_d,1);
    end   
    ftn58temp_d=zeros([ngap_d*16 nw]);
    
    for jj=1:ngap_d
        ftn58temp_d((jj*16-15):(jj*16),2:4)=[
        ipd(jj)+3 ipd(jj)+4      1i*SOd(jj)*2/3            % xy(u)    x2-y2(u)
        ipd(jj)+3 ipd(jj)+2+norb SOd(jj)/3                 % xy(u)    yz(d)
        ipd(jj)+3 ipd(jj)+1+norb -1i*SOd(jj)/3             % xy(u)    xz(d)
        ipd(jj)+2 ipd(jj)+1      1i*SOd(jj)/3              % yz(u)    xz(u)   
        ipd(jj)+2 ipd(jj)+3+norb -SOd(jj)/3                % yz(u)    xy(d)
        ipd(jj)+2 ipd(jj)+4+norb -1i*SOd(jj)/3             % yz(u)    x2-y2(d)
        ipd(jj)+2 ipd(jj)+norb   -1i*SOd(jj)*sqrt(3)/3     % yz(u)    z2(d)
        ipd(jj)+1 ipd(jj)+3+norb 1i*SOd(jj)/3              % xz(u)    xy(d)
        ipd(jj)+1 ipd(jj)+4+norb -SOd(jj)/3                % xz(u)    x2-y2(d)
        ipd(jj)+1 ipd(jj)+norb   sqrt(3)*SOd(jj)/3         % xz(u)    z2(d)
        ipd(jj)+4 ipd(jj)+2+norb 1i*SOd(jj)/3              % x2-y2(u) yz(d)
        ipd(jj)+4 ipd(jj)+1+norb SOd(jj)/3                 % x2-y2(u) xz(d)
        ipd(jj)   ipd(jj)+2+norb 1i*SOd(jj)*sqrt(3)/3      % z2(u)    yz(d)
        ipd(jj)   ipd(jj)+1+norb -SOd(jj)*sqrt(3)/3        % z2(u)    xz(d)
        ipd(jj)+3+norb ipd(jj)+4+norb  -1i*2*SOd(jj)/3     % xy(d)    x2-y2(d)
        ipd(jj)+2+norb ipd(jj)+1+norb  -1i*SOd(jj)/3];     % yz(d)    xz(d)  
    end
    
    ftn58SO         = [zeros([1 nw]); ftn58temp_d];
    ftn58SO(1,1)    = norb*2;
    ftn58SO(1,2)    = ngap_d*16;
    ftn58SO(2:end,1)= (1:ftn58SO(1,2))';
    
%%%--- For both p and d orbitals    
else
    
    ngap_p=length(ipp);
    if length(SOp) == 1
        SOp=repmat(SOp,ngap_p,1);
    end 
    ngap_d=length(ipd);
    if length(SOd) == 1
        SOd=repmat(SOd,ngap_d,1);
    end    
    
    ftn58temp=zeros([ngap_p*6 nw]);
    for ii=1:ngap_p
        ftn58temp((ii*6-5):(ii*6),2:4)= ...
             [ipp(ii) ipp(ii)+1 -1i*SOp(ii)/3            % x(u) y(u)
              ipp(ii) ipp(ii)+norb+2 SOp(ii)/3           % x(u) z(d)
              ipp(ii)+1 ipp(ii)+norb+2 -1i*SOp(ii)/3     % y(u) z(d)
              ipp(ii)+2 ipp(ii)+norb -SOp(ii)/3          % z(u) x(d)
              ipp(ii)+2 ipp(ii)+norb+1 1i*SOp(ii)/3      % z(u) y(d)
              ipp(ii)+norb ipp(ii)+norb+1 1i*SOp(ii)/3]; % x(d) y(d) 
    end

    ftn58temp_d=zeros([ngap_d*16 nw]);
    for jj=1:ngap_d
        ftn58temp_d((jj*16-15):(jj*16),2:4)=[
        ipd(jj)+3 ipd(jj)+4      1i*SOd(jj)*2/3            % xy(u)    x2-y2(u)
        ipd(jj)+3 ipd(jj)+2+norb SOd(jj)/3                 % xy(u)    yz(d)
        ipd(jj)+3 ipd(jj)+1+norb -1i*SOd(jj)/3             % xy(u)    xz(d)
        ipd(jj)+2 ipd(jj)+1      1i*SOd(jj)/3              % yz(u)    xz(u)   
        ipd(jj)+2 ipd(jj)+3+norb -SOd(jj)/3                % yz(u)    xy(d)
        ipd(jj)+2 ipd(jj)+4+norb -1i*SOd(jj)/3             % yz(u)    x2-y2(d)
        ipd(jj)+2 ipd(jj)+norb   -1i*SOd(jj)*sqrt(3)/3     % yz(u)    z2(d)
        ipd(jj)+1 ipd(jj)+3+norb 1i*SOd(jj)/3              % xz(u)    xy(d)
        ipd(jj)+1 ipd(jj)+4+norb -SOd(jj)/3                % xz(u)    x2-y2(d)
        ipd(jj)+1 ipd(jj)+norb   sqrt(3)*SOd(jj)/3         % xz(u)    z2(d)
        ipd(jj)+4 ipd(jj)+2+norb 1i*SOd(jj)/3              % x2-y2(u) yz(d)
        ipd(jj)+4 ipd(jj)+1+norb SOd(jj)/3                 % x2-y2(u) xz(d)
        ipd(jj)   ipd(jj)+2+norb 1i*SOd(jj)*sqrt(3)/3      % z2(u)    yz(d)
        ipd(jj)   ipd(jj)+1+norb -SOd(jj)*sqrt(3)/3        % z2(u)    xz(d)
        ipd(jj)+3+norb ipd(jj)+4+norb  -1i*2*SOd(jj)/3     % xy(d)    x2-y2(d)
        ipd(jj)+2+norb ipd(jj)+1+norb  -1i*SOd(jj)/3];     % yz(d)    xz(d)
    end
    
    ftn58SO         = [zeros([1 nw]); ftn58temp; ftn58temp_d];
    ftn58SO(1,1)    = norb*2;
    ftn58SO(1,2)    = ngap_p*6 + ngap_d*16;
    ftn58SO(2:end,1)= (1:ftn58SO(1,2))';
end

return