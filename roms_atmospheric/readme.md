# 🌤 ROMS Atmospheric Forcing Scripts, continue updating

This repository contains **Python/MATLAB scripts** to download and process atmospheric forcings for the **ROMS** (Regional Ocean Modeling System) model.

---

## 📂 Features

- Automated download of **ERA5 reanalysis** data from the CDS API.
- Support for key atmospheric variables:
  - Wind speed (U10, V10)
  - Air temperature
  - Surface pressure
  - Radiation fluxes
  - Cloud cover
  - Precipitation
  - Relative humidity (from pressure levels)
- Merging and formatting into **ROMS-compatible NetCDF** files.

---

## 📜 Usage
- 3 python scripts (updated 2025/05/17) used to download ERA5: download_ERA5_parallel.py, write_ROMS.py, combine_files.py
  First, determine the date-range/area/variables etc. in "download_ERA5_parallel.py", you may also want to specify the workers in the
  Multipool process. Note: **if you do not want to do parallel download**, just set "with Pool(processes=1) as pool:".

  Then, run it with: python download_ERA5_parallel.py. It will download daily data(e.g.,era5_20010101.zip), unzip and merge the data,
  call a function (write_ROMS.py) to create a ROMS forcing file(e.g., roms_20010101.nc), and put data into the roms file. There will
  be many temporary folders containing the unzipped netcdf files from the downloaded .zip file. You can manually delete all these 
  folders.
  
  After seperate files created, run: python combine_files.py, to merge those files into one forcing file.
  
