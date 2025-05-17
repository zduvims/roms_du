from netCDF4 import Dataset
import numpy as np
import xarray as xr
import os

# ROMS variables
# Uwind: surface u-wind component (m/s)
# Vwind: surface v-wind component (m/s)
# Pair: surface air pressure (mbar)
# Tair: surface air temperature (Kelvin)
# Qair: surface air relative humidity (%)
# rain: rain fall rate (kg/m2/s)
# swrad: net solar shortwave radiation (W/m2)
# lwrad: net olar longwave radiation (W/m2)
# lwrad_down: downward longwave radiation (W/m2)
# cloud: cloud fraction, 0-1, dimensionless




def create_roms_forcing(lon2d, lat2d, time_days, out_filename, output_dir="."):
    """
    Create ROMS forcing NetCDF file with given lon, lat, and time.
    
    lon2d, lat2d: 2D arrays
    time_days: 1D array in Julian days
    out_file: output filename
    """

    # Derive output ROMS filename
    # basename = os.path.basename(era5_filename)  # get file name only
    # date_tag = basename.replace("era5_", "").replace(".nc", "")
    # out_filename = f"roms_{date_tag}.nc"
    out_file = os.path.join(output_dir, out_filename)
    ds_out = Dataset(out_file, "w", format="NETCDF4")

    # Dimensions, the nx should refer to unique longitudes
    # ny refers unique latitudes
    # for example, if the lon_1d has length of 29, lat_1d has length of 33, then
    # lon2d/lat2d should have size of 33*29
    ny,nx = lon2d.shape
    nt = len(time_days)
    ds_out.createDimension("xrho", nx)
    ds_out.createDimension("yrho", ny)
    ds_out.createDimension("time", nt)

    # Coordinates
    time_var = ds_out.createVariable("time", "f8", ("time",))
    time_var.long_name = "atmospheric forcing time"
    time_var.units = "days since 1858-11-17"
    time_var.field = "time, scalar, series"
    time_var[:] = time_days

    lon_var = ds_out.createVariable("lon", "f8", ("yrho", "xrho"))
    lon_var.long_name = "longitude"
    lon_var.units = "degrees_east"
    lon_var.field = "xp, scalar, series"
    lon_var[:, :] = lon2d

    lat_var = ds_out.createVariable("lat", "f8", ("yrho", "xrho"))
    lat_var.long_name = "latitude"
    lat_var.units = "degrees_north"
    lat_var.field = "yp, scalar, series"
    lat_var[:, :] = lat2d

    # === Add fixed variables ===
    # make sure you check "varinfo.dat" in the ROMS source code, the defined var name and unit must
    # be the same as in varinfo.dat
    
    # Example for Uwind/Vwind:
    ds_out.createDimension("wind_time", nt)
    wind_time = ds_out.createVariable("wind_time", "f8", ("wind_time"))
    wind_time.units = "days since 1858-11-17"
    wind_time[:] = time_days

    Uwind = ds_out.createVariable("Uwind", "f8", ("wind_time", "yrho", "xrho" ))
    Uwind.units = "m/s"
    Uwind.long_name = "surface u-wind component"
    Uwind.coordinates = "lon lat"
    Uwind.time = "wind_time"
    
    Vwind = ds_out.createVariable("Vwind", "f8", ("wind_time", "yrho", "xrho" ))
    Vwind.units = "m/s"
    Vwind.long_name = "surface v-wind component"
    Vwind.coordinates = "lon lat"
    Vwind.time = "wind_time"

    # cloud fraction, 0-1, dimensionless
    ds_out.createDimension("cloud_time", nt)
    cloud_time = ds_out.createVariable("cloud_time", "f8", ("cloud_time"))
    cloud_time.units = "days since 1858-11-17"
    cloud_time[:] = time_days
    cloud = ds_out.createVariable("cloud", "f8", ("cloud_time","yrho", "xrho" ))
    cloud.units = "nondimensional"
    cloud.long_name = "cloud fraction"
    cloud.coordinates = "lon lat"
    cloud.time = "cloud_time"

    # air temperature, Tair, Kelvin
    ds_out.createDimension("Tair_time", nt)
    Tair_time = ds_out.createVariable("Tair_time", "f8", ("Tair_time"))
    Tair_time.units = "days since 1858-11-17"
    Tair_time[:] = time_days
    Tair = ds_out.createVariable("Tair", "f8", ("Tair_time", "yrho", "xrho" ))
    Tair.units = "kelvin"
    Tair.long_name = "2m surface air temperature"
    Tair.coordinates = "lon lat"
    Tair.time = "Tair_time"

    # net surface long-wave radiation flux, W/m2
    ds_out.createDimension("lwrad_time", nt)
    lwrad_time = ds_out.createVariable("lwrad_time", "f8", ("lwrad_time"))
    lwrad_time.units = "days since 1858-11-17"
    lwrad_time[:] = time_days
    lwrad = ds_out.createVariable("lwrad", "f8", ("lwrad_time", "yrho", "xrho" ))
    lwrad.units = "W/m2"
    lwrad.long_name = "mean_surface_net_long_wave_radiation_flux"
    lwrad.coordinates = "lon lat"
    lwrad.time = "lwrad_time"

    # downwelling surface long-wave radiation flux, W/m2
    lwrad_down = ds_out.createVariable("lwrad_down", "f8", ("lwrad_time", "yrho", "xrho"))
    lwrad_down.units = "W/m2"
    lwrad_down.long_name = "mean_surface_downward_long_wave_radiation_flux"
    lwrad_down.coordinates = "lon lat"
    lwrad_down.time = "lwrad_time"

    # net surface short-wave radiation flux, W/m2
    ds_out.createDimension("swrad_time", nt)
    swrad_time = ds_out.createVariable("swrad_time", "f8", ("swrad_time"))
    swrad_time.units = "days since 1858-11-17"
    swrad_time[:] = time_days
    swrad = ds_out.createVariable("swrad", "f8", ( "swrad_time", "yrho", "xrho"))
    swrad.units = "W/m2"
    swrad.long_name = "mean_surface_net_short_wave_radiation_flux"
    swrad.coordinates = "lon lat"
    swrad.time = "swrad_time"

    # precipitation rate, kg/m2/s
    ds_out.createDimension("rain_time", nt)
    rain_time = ds_out.createVariable("rain_time", "f8", ("rain_time"))
    rain_time.units = "days since 1858-11-17"
    rain_time[:] = time_days
    rain = ds_out.createVariable("rain", "f8", ("rain_time", "yrho", "xrho" ))
    rain.units = "kg/m2/s"
    rain.long_name = "Time-mean total precipitation rate"
    rain.coordinates = "lon lat"
    rain.time = "rain_time"

    # surface air pressure, millibar, Pair
    ds_out.createDimension("Pair_time", nt)
    Pair_time = ds_out.createVariable("Pair_time", "f8", ("Pair_time"))
    Pair_time.units = "days since 1858-11-17"
    Pair_time[:] = time_days
    Pair = ds_out.createVariable("Pair", "f8", ("Pair_time", "yrho", "xrho" ))
    Pair.units = "millibar"
    Pair.long_name = "air_pressure_at_mean_sea_level"
    Pair.coordinates = "lon lat"
    Pair.time = "Pair_time"

    # relative humidity, Qair, percentage
    ds_out.createDimension("Qair_time", nt)
    Qair_time = ds_out.createVariable("Qair_time", "f8", ("Qair_time"))
    Qair_time.units = "days since 1858-11-17"
    Qair_time[:] = time_days
    Qair = ds_out.createVariable("Qair", "f8", ("Qair_time","yrho", "xrho" ))
    Qair.units = "percent"
    Qair.long_name = "surface air relative humidity"
    Qair.coordinates = "lon lat"
    Qair.time = "Qair_time"

    ds_out.close()
    return out_file


