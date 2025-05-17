import cdsapi
import zipfile
import xarray as xr
from netCDF4 import Dataset, num2date
import os
from roms_atmospheric.write_ROMS import create_roms_forcing
from datetime import datetime, timedelta
import numpy as np
from multiprocessing import Pool
import time


# Output directory
output_dir = "/sciclone/data10/zdu/projects/Andaman/files/forcings/ERA5/year2005/"
os.makedirs(output_dir, exist_ok=True)

# Log file for failed downloads
log_file = os.path.join(output_dir, "download_failed_log.txt")

# Date range you want to grab the data
start_date = datetime(2005, 7, 1)
end_date   = datetime(2005, 12, 31)

# Define your variables and region
variables = [
    "10m_u_component_of_wind",
    "10m_v_component_of_wind",
    "2m_dewpoint_temperature",
    "2m_temperature",
    "mean_sea_level_pressure",
    "mean_surface_downward_long_wave_radiation_flux",
    "mean_surface_net_long_wave_radiation_flux",
    "mean_surface_net_short_wave_radiation_flux",
    "mean_total_precipitation_rate",
    "total_cloud_cover"
]
area = [20, 92, 12, 99]  # North, West, South, East

# Dataset to download
dataset = "reanalysis-era5-single-levels"

# loop download for the date range 
# files downloaded from ERA5 is a zip format, which contains multiple nc files
# This routine will download the ERA5 data, unzip it, combines variables from the nc files and
# save it to a new nc file. Finally, create a ROMS forcing nc file and put data into it.
dates = []
d = start_date
while d <= end_date:
    dates.append(d)
    d += timedelta(days=1)

def download_and_process(date):
    
    # Initialize CDS API client
    client = cdsapi.Client()

    year_str = str(date.year)
    month_str = f"{date.month:02d}"
    day_str = f"{date.day:02d}"
    date_tag = f"{year_str}{month_str}{day_str}"

    # define the downloaded zip filename
    zip_file = os.path.join(output_dir, f"era5_{year_str}{month_str}{day_str}.zip")
    # define the combined ERA5 nc file
    combined_file = os.path.join(output_dir, f"era5_{year_str}{month_str}{day_str}.nc")
    # define temporary directory, because each zip file has the same nc filenames, if
    # multiple workers unzip at the same time, it will cause serious problem
    tmp_dir = os.path.join(output_dir, f"tmp_{date_tag}")
    # Skip if already done
    if os.path.exists(zip_file) or os.path.exists(combined_file) or os.path.exists(tmp_dir):
          print(f"✅ Already exists: {zip_file}")
          return
    
    os.makedirs(tmp_dir, exist_ok=True)
    
    # EAR5 request
    request = {
        "product_type": "reanalysis",
        "variable": variables,
        "year": year_str,
        "month": month_str,
        "day": day_str,
        "time": [f"{h:02d}:00" for h in range(24)],
        "area": area,
        "format": "netcdf"
    }

    try:
        print(f"⬇️ Downloading {zip_file}")
        client.retrieve(dataset, request).download(zip_file)
        print(f"✅ Finished: {zip_file}")

        # unzip and combine, then delete unzipped nc files
        print(f"Unzipping {zip_file}")
        # for each zip_file, create a temporary directort and zip the downloaded file
        with zipfile.ZipFile(zip_file,'r') as zip_ref:
             zip_ref.extractall(tmp_dir)
        
        # combine the unzipped files
        # Merge extracted files
        file_list = [os.path.join(tmp_dir, f) for f in zip_ref.namelist()]
        ds_merged = xr.merge([xr.open_dataset(f) for f in file_list])
        ds_merged.to_netcdf(combined_file)
        
        #%% Download relative humidity from pressure level dataset (at 1000 hPa)
        rh_nc_file = os.path.join(output_dir, f"era5_rh_{year_str}{month_str}{day_str}.nc")

        if not os.path.exists(rh_nc_file):
            try:
                rh_request = {
                    "product_type": "reanalysis",
                    "variable": ["relative_humidity"],
                    "pressure_level": ["1000"],
                    "year": year_str,
                    "month": month_str,
                    "day": day_str,
                    "time": [f"{h:02d}:00" for h in range(24)],
                    "area": area,
                    "format": "netcdf"
                }
                print(f"⬇️ Downloading RH: {rh_nc_file}")
                client.retrieve("reanalysis-era5-pressure-levels", rh_request).download(rh_nc_file)
                print(f"✅ Finished RH: {rh_nc_file}")

            except Exception as e:
                print(f"❌ Failed RH: {rh_nc_file} — {e}")
                with open(os.path.join(output_dir, log_file), 'a') as f:
                    f.write(f"{year_str}-{month_str}-{day_str} RH: {str(e)}\n")
            
    # if donwload failed
    except Exception as e:
        print(f"❌ Failed: {zip_file} — {e}")
        with open(os.path.join(output_dir, log_file), 'a') as f:
             f.write(f"{year_str}-{month_str}-{day_str}: {str(e)}\n")
    
    #%% write to ROMS
    ds = Dataset(combined_file)
    lon = ds.variables['longitude'][:]
    lat = ds.variables['latitude'][:]
    lon2d,lat2d = np.meshgrid(lon,lat, indexing="xy")
    time = ds.variables['valid_time'] # seconds since 1970-1-1
    # convert to julian days
    datetime_vals = num2date(time[:], units=time.units)
    julian_base = datetime(1858,11,17)
    julian_days = [(dt - julian_base).total_seconds() / 86400 for dt in datetime_vals]

    # call to create the ROMS forcing file
    # Derive output ROMS filename
    basename = os.path.basename(combined_file)  # get file name only
    date_tag = basename.replace("era5_", "").replace(".nc", "")
    roms_filename = f"roms_{date_tag}.nc"
    out_file = os.path.join(output_dir, roms_filename)
    if not (
        os.path.exists(combined_file) and
        os.path.exists(rh_nc_file) and
        not os.path.exists(out_file)
    ):
        return  # Skip if any condition fails
        
    create_roms_forcing(lon2d, lat2d, julian_days, roms_filename, output_dir)

    # Re-open ROMS forcing file to insert data
    # ROMS variables has dimension of (time,lat,lon)
    # ds variables has dimension of (time,lat,lon)
    with Dataset(os.path.join(output_dir, roms_filename), 'a') as roms_ds:
         roms_ds['Uwind'][:,:,:] = ds.variables['u10'][:]
         roms_ds['Vwind'][:,:,:] =ds.variables['v10'][:]

         # surface temperature
         roms_ds['Tair'][:,:,:] = ds.variables['t2m'][:]
         # mean surface pressure, ROMS needs millibar, while data is Pa
         roms_ds['Pair'][:,:,:] = ds.variables['msl'][:] / 100.0
         # cloud fraction
         roms_ds['cloud'][:,:,:] = ds.variables['tcc'][:]
         # downward long-wave radiation flux
         roms_ds['lwrad_down'][:,:,:] = ds.variables['avg_sdlwrf'][:]
         # net long-wave radiation flux
         roms_ds['lwrad'][:,:,:] = ds.variables['avg_snlwrf'][:]
         # net short-wave radiation flux
         roms_ds['swrad'][:,:,:] = ds.variables['avg_snswrf'][:]
         # rainfall rate
         roms_ds['rain'][:,:,:] = ds.variables['avg_tprate'][:]
         # relative humidity, it has dimension of (time,1,lat,lon)
         rh_ds_nc = Dataset(rh_nc_file)
         rh_raw = rh_ds_nc.variables['r'][:]
         #rh_squeezed = np.squeeze(rh_raw, axis=1)
         roms_ds['Qair'][:, :, :] = np.squeeze(rh_raw, axis=1)
         rh_ds_nc.close()
        
    # remove ERA5 nc file after ROMS file created
    os.remove(zip_file)
    os.remove(combined_file)
    os.remove(rh_nc_file)


# Use multiprocessing Pool
if __name__ == "__main__":

    start_time = time.time()
    with Pool(processes=10) as pool:  # set number of workers here
        pool.map(download_and_process, dates)
    
    end_time = time.time()
    total_seconds = end_time - start_time

    print("✅ All files downloaded.")
    print(f"🕒 Total runtime: {total_seconds:.2f} seconds")





