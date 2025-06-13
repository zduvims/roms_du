function [fn]=updatclim_coawst_mw_parallel_zdu(T1, gn, clm, clmname, wdr, url)
% Modified by Brandy Armstrong January 2012 to use only NCTOOLBOX 
% and Matlab builtin functions to read and write netcdf files
% jcw Feb 2019 - only use matalb BI
%
%T1 = date for climatology file
%gn = data from grid
%clm = data of hycom indices
%wdr = the working directory
%clmname = grid name prefix for climatology filenames
%url = where get data from

%
%determine indices for time period of interpolation
%
disp('getting the number of time records ...');
try 
  time=ncread(url,'MT');
  t0=datenum(1900,12,31); % tr0=datenum(1858,11,17);
  tg=time+t0;
catch
  time=ncread(url,'time')/24;
  time_att=ncreadatt(url,'time','units'); %'hours since 2000-01-01 00:00:00'
  tg=datenum(str2num(time_att(13:16)),str2num(time_att(18:19)),str2num(time_att(21:22)),str2num(time_att(24:25)),str2num(time_att(27:28)),str2num(time_att(30:31)));
  tg=time+tg;
end
tg2=julian(str2num(datestr(tg,'yyyy')),str2num(datestr(tg,'mm')),str2num(datestr(tg,'dd')),str2num(datestr(tg,'HH')))-2400001;

% get user times
t1 = tg - T1;
f = find(t1 >= 0 & t1 < 1);
if isempty(f)
  error(['HYCOM file does not cover the date',...
      ' at which you want to create the boundary file', newline,...
      datestr(T1)]);
end
tid1 = f(1); 

fn=clmname;
disp(['creating netcdf file ',fn]);
create_roms_netcdf_clm_mwUL(fn,gn,1);% converted to BI functions

%fill grid dims using builtin (BI) functions
RN=netcdf.open(fn,'NC_WRITE');
lonid=netcdf.inqVarID(RN,'lon_rho');
netcdf.putVar(RN,lonid,gn.lon_rho);
latid=netcdf.inqVarID(RN,'lat_rho');
netcdf.putVar(RN,latid,gn.lat_rho);

tempid=netcdf.inqVarID(RN,'ocean_time');
netcdf.putVar(RN,tempid,tg2(tid1));
tempid=netcdf.inqVarID(RN,'zeta_time');
netcdf.putVar(RN,tempid,tg2(tid1));
tempid=netcdf.inqVarID(RN,'v2d_time');
netcdf.putVar(RN,tempid,tg2(tid1));
tempid=netcdf.inqVarID(RN,'v3d_time');
netcdf.putVar(RN,tempid,tg2(tid1));
tempid=netcdf.inqVarID(RN,'salt_time');
netcdf.putVar(RN,tempid,tg2(tid1));
tempid=netcdf.inqVarID(RN,'temp_time');
netcdf.putVar(RN,tempid,tg2(tid1));
netcdf.close(RN)

%% interpolate 3D u/v from HYCOM to ROMS grid
tz_levs=length(clm.z);
X=clm.lon;
Y=clm.lat;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Interpolating u for ',datestr(tg(tid1))]);
ttu=1;
clm_u=zeros([size(gn.lon_rho) length(clm.z)]); % modify by zdu
while ttu==1
    try
        try
          tmpt=ncread(url,'u',[clm.ig0 clm.jg0 1 tid1],...
              [clm.ig1-clm.ig0+1 clm.jg1-clm.jg0+1 tz_levs 1 ] );
        catch
          tmpt=ncread(url,'water_u',[clm.ig0 clm.jg0 1 tid1],...
              [clm.ig1-clm.ig0+1 clm.jg1-clm.jg0+1 tz_levs 1 ] );
        end
        parfor k=1:tz_levs
            disp(['doing griddata u for HYCOM level ' num2str(k)]);
            tmp=double(squeeze(tmpt(:,:,k)));
            if (k==1)
              Fu = scatteredInterpolant(X(:),Y(:),tmp(:));
            else
              Fu.Values = tmp(:);
            end
            cff = Fu(gn.lon_rho,gn.lat_rho);
            clm_u(:,:,k)=maplev(cff);
        end
        ttu=0;
    catch
        disp(['catch u Unable to download HYCOM u data at' datestr(tg(tid1))]);
        fid=fopen('coawstlog.txt','a');
        fprintf(fid,'Unable to download HYCOM u data at ');
        fprintf(fid,datestr(tg(tid1)));
        fprintf(fid,'\n');
    end
end
%== Vertical interpolation (t,s,u,v) from standard z-level to s-level
u=roms_from_stdlev_mw_zdu(gn.lon_rho,gn.lat_rho,clm.z,clm_u,gn,'u',0); %% modify by zdu
%clm=rmfield(clm,'u');
save u.mat u
%clear u;

disp(['Interpolating v for ',datestr(tg(tid1))]);
ttv=1;
clm_v=zeros([size(gn.lon_rho) length(clm.z)]);
while ttv==1
    try
        try
          tmpt=ncread(url,'v',[clm.ig0 clm.jg0 1 tid1],...
              [clm.ig1-clm.ig0+1 clm.jg1-clm.jg0+1 tz_levs 1 ] );
        catch
          tmpt=ncread(url,'water_v',[clm.ig0 clm.jg0 1 tid1],...
              [clm.ig1-clm.ig0+1 clm.jg1-clm.jg0+1 tz_levs 1 ] );
        end
        parfor k=1:tz_levs
            disp(['doing griddata v for HYCOM level ' num2str(k)]);
            tmp=double(squeeze(tmpt(:,:,k)));
            if (k==1)
              Fv = scatteredInterpolant(X(:),Y(:),tmp(:));
            else
              Fv.Values = tmp(:);
            end
            cff = Fv(gn.lon_rho,gn.lat_rho);
            clm_v(:,:,k)=maplev(cff); % modify by zdu
        end
        ttv=0;
    catch
        disp(['catch v Unable to download HYCOM v data at' datestr(tg(tid1))]);
        fid=fopen('coawstlog.txt','a');
        fprintf(fid,'Unable to download HYCOM v data at');
        fprintf(fid,datestr(tg(tid1)));
        fprintf(fid,'\n');
    end
end
%== Vertical interpolation (t,s,u,v) from standard z-level to s-level
v=roms_from_stdlev_mw_zdu(gn.lon_rho,gn.lat_rho,clm.z,clm_v,gn,'v',0); %% modify by zdu
%clm=rmfield(clm,'v');
save v.mat v
%clear v;

%== Rotate the velocity
theta=exp(-sqrt(-1)*mean(mean(gn.angle)));
%load u.mat; load v.mat
disp('doing rotation to grid for u and v');
uv=(u2rho_3d_mw(u)+sqrt(-1)*v2rho_3d_mw(v)).*theta;
u=rho2u_3d_mw(real(uv)); v=rho2v_3d_mw(imag(uv));
clear uv

% write outputs to the climatology file
RN=netcdf.open(fn,'NC_WRITE');
tempid=netcdf.inqVarID(RN,'u');
netcdf.putVar(RN,tempid,u);
tempid=netcdf.inqVarID(RN,'v');
netcdf.putVar(RN,tempid,v);
netcdf.close(RN);

%% depth average 3D u,v to get ubar/vbar
% u/v already rotated, no need again here
cc=roms_zint_mw_zdu(u,gn);  ubar=rho2u_2d_mw(u2rho_2d_mw(cc)./gn.h);
cc=roms_zint_mw_zdu(v,gn);  vbar=rho2v_2d_mw(v2rho_2d_mw(cc)./gn.h);

RN=netcdf.open(fn,'NC_WRITE');
tempid=netcdf.inqVarID(RN,'ubar');
netcdf.putVar(RN,tempid,ubar);
tempid=netcdf.inqVarID(RN,'vbar');
netcdf.putVar(RN,tempid,vbar);
netcdf.close(RN);

clear ubar u clm_u
clear vbar v clm_v

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% interpolate the zeta data
disp(['Interpolating zeta for ',datestr(tg(tid1))]);
ttz=1;
while ttz==1
    try
        try
          tmpt=ncread(url,'ssh',[clm.ig0 clm.jg0 tid1],...
              [clm.ig1-clm.ig0+1 clm.jg1-clm.jg0+1 1 ] );
        catch
          tmpt=ncread(url,'surf_el',[clm.ig0 clm.jg0 tid1],...
              [clm.ig1-clm.ig0+1 clm.jg1-clm.jg0+1 1 ] );
        end
        tmp=double(squeeze(tmpt(:,:)));
        disp(['doing griddata zeta for HYCOM at ', datestr(tg(tid1))]);
        F = scatteredInterpolant(X(:),Y(:),tmp(:));
        cff = F(gn.lon_rho,gn.lat_rho);
        zeta=maplev(cff);
        ttz=0;
    catch
        disp(['catch z Unable to download HYCOM ssh data at' datestr(tg(tid1))]);
        fid=fopen('coawstlog.txt','a');
        fprintf(fid,'Unable to download HYCOM ssh data at');
        fprintf(fid,datestr(tg(tid1)));
        fprintf(fid,'\n');
    end
end
clear tmp

% write zeta to climatology file
RN=netcdf.open(fn,'NC_WRITE');
tempid=netcdf.inqVarID(RN,'zeta');
netcdf.putVar(RN,tempid,zeta);
netcdf.close(RN);
clear zeta


%% Interpolate temperature from HYCOM to ROMS grid
disp(['Interpolating temp for ',datestr(tg(tid1))]);
ttt=1;
clm_temp=zeros([size(gn.lon_rho) length(clm.z)]);
while ttt==1
    try
        try
          tmpt=ncread(url,'temperature',[clm.ig0 clm.jg0 1 tid1],...
              [clm.ig1-clm.ig0+1 clm.jg1-clm.jg0+1 tz_levs 1 ] );
        catch
          tmpt=ncread(url,'water_temp',[clm.ig0 clm.jg0 1 tid1],...
              [clm.ig1-clm.ig0+1 clm.jg1-clm.jg0+1 tz_levs 1 ] );
        end

        parfor k=1:tz_levs
            disp(['doing griddata temp for HYCOM level ' num2str(k)]);
            tmp=double(squeeze(tmpt(:,:,k)));
            if (k==1)
              Ft = scatteredInterpolant(X(:),Y(:),tmp(:));
            else
              Ft.Values = tmp(:);
            end
            cff = Ft(gn.lon_rho,gn.lat_rho);
%           cff(cff<0)=nan;
            clm_temp(:,:,k)=maplev(cff);
        end
        ttt=0;
    catch
        disp(['catch temp Unable to download HYCOM temp data at' datestr(tg(tid1))]);
        fid=fopen('coawstlog.txt','a');
        fprintf(fid,'Unable to download HYCOM temp data at ');
        fprintf(fid,datestr(tg(tid1)));
        fprintf(fid,'\n');
    end
end

%== Vertical interpolation (t,s,u,v) from standard z-level to s-level
temp=roms_from_stdlev_mw_zdu(gn.lon_rho,gn.lat_rho,clm.z,clm_temp,gn,'rho',0);
% clm=rmfield(clm,'temp');

% write temp to climatology file
RN=netcdf.open(fn,'NC_WRITE');
tempid=netcdf.inqVarID(RN,'temp');
netcdf.putVar(RN,tempid,temp);
netcdf.close(RN);
clear temp clm_temp


%% Interpolate salinity from HYCOM to ROMS grid
disp(['Interpolating salt for ',datestr(tg(tid1))]);
tts=1;
clm_salt=zeros([size(gn.lon_rho) length(clm.z)]);
while tts==1
    try
        tmpt=ncread(url,'salinity',[clm.ig0 clm.jg0 1 tid1],...
            [clm.ig1-clm.ig0+1 clm.jg1-clm.jg0+1 tz_levs 1 ] );

        parfor k=1:tz_levs
            disp(['doing griddata salt for HYCOM level ' num2str(k)]);
            tmp=double(squeeze(tmpt(:,:,k)));
            if (k==1)
              Fs = scatteredInterpolant(X(:),Y(:),tmp(:));
            else
              Fs.Values = tmp(:);
            end
            cff = Fs(gn.lon_rho,gn.lat_rho);
            cff(cff<0)=nan;
            clm_salt(:,:,k)=maplev(cff);
        end
        tts=0;
    catch
        disp(['catch temp Unable to download HYCOM temp data at' datestr(tg(tid1))]);
        fid=fopen('coawstlog.txt','a');
        fprintf(fid,'Unable to download HYCOM temp data at ');
        fprintf(fid,datestr(tg(tid1)));
        fprintf(fid,'\n');
    end
end

%== Vertical interpolation (t,s,u,v) from standard z-level to s-level
salt=roms_from_stdlev_mw_zdu(gn.lon_rho,gn.lat_rho,clm.z,clm_salt,gn,'rho',0);
%clm=rmfield(clm,'salt');

% write salt to climatology file
RN=netcdf.open(fn,'NC_WRITE');
tempid=netcdf.inqVarID(RN,'salt');
netcdf.putVar(RN,tempid,salt);
netcdf.close(RN);
clear salt clm_salt

disp(['Finished creating clim file at ' datestr(tg(tid1))]);
%%
