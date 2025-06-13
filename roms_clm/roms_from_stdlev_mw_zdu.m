function roms = roms_from_stdlev_mw_zdu(lon,lat,zlev,data,grd,CgridPos,do_interp2)
% Interpolate a 3D gridded (e.g. climatology or synthetic t(z),s(z) analysis)
% to a 3D ROMS grid.
% 
% Input "data" must be dimensioned (x, y, z, time)
%
% roms = roms_from_stdlev_mw(lon,lat,zlev,data,grd,CgridPos,do_interp2)
%
% Zhiyun Du, modified version (2025)
%   -- many changes have been made, the original "roms_from_stdlev_mw" assumes the input data
%   -- has the dimension (time,x,y,z), I changed the code to read data in the dimension of (x, y, z, time)

if nargin < 6
  CgridPos = 'rho';
end

switch CgridPos
  case 'rho'
    romslon = grd.lon_rho;
    romslat = grd.lat_rho;
    romsmask = grd.mask_rho;
  case 'u'
    romslon = grd.lon_u;
    romslat = grd.lat_u;
    romsmask = grd.mask_u;
  case 'v'
    romslon = grd.lon_v;
    romslat = grd.lat_v;
    romsmask = grd.mask_v;
  otherwise
    error('Input CgridPos must be ''rho'', ''u'' or ''v'' ')
end

land = find(romsmask == 0);

roms = [];

if isempty(zlev)
  if ndims(data) == 4
    Nt = size(data,4);
    roms_stdlev = zeros([size(romslon), Nt]);
  elseif ndims(data) == 3
    Nt = 1;
    roms_stdlev = zeros([size(romslon)]);
  else
    error('Check the size of "data" ')
  end
else
  if ndims(data) == 4
    Nt = size(data,4);
    Nz = size(data,3);
    roms_stdlev = zeros([size(romslon) Nz Nt]);
    Ns = size(grd.z_r,3); %% make sure the z variables has the size (x,y,z)
    roms = zeros([size(romslon) Nz Nt]);
  elseif ndims(data) == 3
    Nt = 1;
    Nz = size(data,3);
    roms_stdlev = zeros([size(romslon) Nz]);
    Ns = size(grd.z_r,3);
    roms = zeros([size(romslon), Ns]);
  elseif ndims(data) == 2
    Nt = 1;
    Nz = 1;
    roms_stdlev = zeros([size(romslon)]);
  else
    error('Check the size of "data" ')
  end

  % force zlev to be negative values to be consistent with ROMS
  % convention (i.e. trap the case of z given as standard depths (>0) 
  zlev = -abs(zlev(:));

  if Nz > 1
    % check that the data is ordered from deep to shallow
    if any(diff(zlev) < 0)
        % depths are arranged from shallowest to deepest, so flip
      warning('Reversing zlev to arrange as deep to shallow');
      zlev = flip(zlev);
      data = flip(data,3);
    end
  end
end


% INTERPOLATE TO HORIZONTAL GRID ----------------------------------------------
% -----------------------------------------------------------------------------

if nargin < 7
  do_interp2 = 1;
end


if do_interp2
  disp('Interpolating to the ROMS horizontal grid at standard levels')
  switch ndims(roms_stdlev)
    case 4
      for l = 1:Nt
        disp([' Doing time ' int2str(l) ' of ' int2str(Nt)])
        for k = 1:Nz
          datawrk = squeeze(data(:,:,k,l));
          tmp = interp2(lon,lat,datawrk,romslon,romslat,'spline');
          if ~isempty(land), tmp(land) = 0; end
          roms_stdlev(:,:,k,l) = tmp;
        end
      end
    case 3
      for k = 1:max([Nt Nz])
        datawrk = squeeze(data(:,:,k));
        try
          tmp = interp2(lon,lat,datawrk,romslon,romslat,'spline');
          if ~isempty(land), tmp(land) = 0; end
          roms_stdlev(:,:,k) = tmp;
        catch
          tmp = griddata(lon,lat,datawrk,romslon,romslat);
          if ~isempty(land), tmp(land) = 0; end
          roms_stdlev(:,:,k) = tmp;
        end
      end
    case 2
      datawrk = data;
      try
        tmp = interp2(lon,lat,datawrk,romslon,romslat,'spline');
        if ~isempty(land), tmp(land) = 0; end
        roms_stdlev(:,:) = tmp;
      catch
        tmp = griddata(lon,lat,datawrk,romslon,romslat);
        if ~isempty(land), tmp(land) = 0; end
        roms_stdlev(:,:) = tmp;
      end
  end

else

  disp('Input is assumed to be on the ROMS horizontal rho grid')
  switch CgridPos
    case 'rho'
      roms_stdlev = data;
    case 'u'
      disp(' but is being averaged to the u-points grid ')
      [LL,MM] = size(romslon);
      switch ndims(data)
        case 4
          roms_stdlev = 0.5*(data(1:LL,:,:,:) + data(2:LL+1,:,:,:));
        case 3
          roms_stdlev = 0.5*(data(1:LL,:,:) + data(2:LL+1,:,:));
        case 2
          roms_stdlev = 0.5*(data(1:LL,:) + data(2:LL+1,:));
      end

    case 'v'
      disp(' but is being averaged to the v-points grid ')
      [LL,MM] = size(romslon);
      switch ndims(data)
        case 4
          roms_stdlev = 0.5*(data(:,1:MM,:,:) + data(:,2:MM+1,:,:));
        case 3
          roms_stdlev = 0.5*(data(:,1:MM,:) + data(:,2:MM+1,:));
        case 2
          roms_stdlev = 0.5*(data(:,1:MM) + data(:,2:MM));
      end

  end
end

if isempty(roms)
  disp('No vertical interpolation is required')
  roms = roms_stdlev;
else
  disp('Interpolating to ROMS s-coordinates')

  % if necessary, average the roms z_r to velocity points
  z_ = grd.z_r;
  switch CgridPos
    case 'u'
      L = size(z_,1);
      z_ = 0.5*(z_(1:L-1,:,:) + z_(2:L,:,:));
    case 'v'
      M = size(z_,2);
      z_ = 0.5*(z_(:,1:M-1,:) + z_(:,2:M,:));
  end

  switch ndims(roms)
    case 4
      for l = 1:Nt
        disp([' Doing time ' int2str(l) ' of ' int2str(Nt)])
        % interpolate a y-z plane each time
        Nx = size(roms,1);
        for i = 1:Nx
          if ~rem(i,20)
            disp(['  Doing i = ' int2str(i) ' of ' int2str(Nx)])
          end
          z = squeeze(z_(i,:,:));
          M = size(z,1);
          x = repmat((1:M)',[1 Ns]);

          % There may be ROMS z values outside the stdlev z range, so pad
          % above and below before interp2 (just like in roms_zslice)
          % (Hmmm ... there may still be a catch if there are some very
          % deep depths with NaNs in the data)
          [za,xa] = meshgrid([-10000; -abs(zlev); 10],1:M);
          data_interp = squeeze(roms_stdlev(i,:,:,l));
          data_interp = [data_interp(:,1), data_interp, data_interp(:,end)];
          roms(i,:,:,l) = interp2(xa',za',data_interp',x,z,'spline');
        end
      end
    case 3
      clear lat lon romslat romslon romsmask
      Nx = size(roms,1);
      for i = 1:Nx
        if ~rem(i,20)
          disp(['  Doing i = ' int2str(i) ' of ' int2str(Nx)])
        end
        z = squeeze(z_(i,:,:));
        M = size(z,1);
        x = repmat((1:M)',[1 Ns]);
        [za,xa] = meshgrid([-10000; -abs(zlev); 10],1:M);
        data_interp = squeeze(roms_stdlev(i,:,:));
        data_interp = [data_interp(:,1), data_interp, data_interp(:,end)];
        data_interp((isnan(data_interp)==1)) = 0;
        roms(i,:,:) = interp2(xa',za',data_interp',x,z,'spline');
      end
  end
end
end
