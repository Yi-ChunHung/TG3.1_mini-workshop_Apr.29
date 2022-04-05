function lat = latticegen(BR,nx,ny,nz)
%% Generate the positions of cells in Cartesian coordinates 

    [x, y, z] = ndgrid(-nx:nx, -ny:ny, -nz:nz);
    xyz       = [x(:), y(:), z(:)];
    lat       = (xyz*BR)';

%lat=[];
%for ii=-nx:nx
%    x=ii*BR(:,1);
%    for jj=-ny:ny
%        y=jj*BR(:,2);
%        for kk=-nz:nz
%            lat=[lat x+y+kk*BR(:,3)];
%            % lat = [lat BR*[ii;jj;kk]];
%        end
%    end
%end
