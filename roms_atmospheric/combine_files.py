import xarray as xr
import os
from netCDF4 import Dataset
import numpy as np
from write_ROMS import create_roms_forcing

# Path to your downloaded seperate nc files
input_dir = "/sciclone/data10/zdu/projects/Andaman/files/forcings/ERA5/year2007/"  # change this
output_file = "roms_ERA5_2007.nc"    # output filename
output_path = os.path.join(input_dir, output_file)

# Get sorted list of files (very important to preserve chronological order)
roms_files = sorted([
    os.path.join(input_dir, f)
    for f in os.listdir(input_dir)
    if f.startswith("roms") and f.endswith(".nc")
])

# Read lat/lon from the first file
with Dataset(roms_files[0]) as ds:
    lon2d = ds.variables["lon"][:, :]
    lat2d = ds.variables["lat"][:, :]

# get total time records
# Step 2: Gather actual "time" values
all_times = []
for f in roms_files:
    with Dataset(f) as ds:
        all_times.append(ds.variables["time"][:])
all_times = np.concatenate(all_times)

# Create combined output file
create_roms_forcing(lon2d, lat2d, all_times, output_file, output_dir=input_dir)

# Open the combined output file
variables = [
    "Uwind", "Vwind", "Tair", "Pair", "Qair", "cloud",
    "lwrad", "lwrad_down", "swrad", "rain"
]
time_vars = {
    "Uwind": "wind_time", "Vwind": "wind_time",
    "Tair": "Tair_time", "Pair": "Pair_time", "Qair": "Qair_time",
    "cloud": "cloud_time", "lwrad": "lwrad_time",
    "lwrad_down": "lwrad_time", "swrad": "swrad_time",
    "rain": "rain_time"
}

with Dataset(output_path, "a") as ds_out:
    time_index = 0
    for f in roms_files:
        with Dataset(f) as ds_in:
            nt = ds_in.dimensions["time"].size
            for var in variables:
                tvar = time_vars[var]
                ds_out[var][time_index:time_index+nt, :, :] = ds_in[var][:]

            time_index += nt

print(f"✅ All ROMS forcing variables filled into {output_file}")
