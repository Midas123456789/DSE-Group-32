import numpy as np
import xarray as xr
import argparse
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from tqdm import tqdm
import matplotlib.colors # For setting NaN color

# File: b00885fac425ddae6f2421ac6486807a.nc
# Product Type: Reanalysis
# Variables: U-component of wind, V-component of wind
# Year of Data: 2022, 2023, 2024
# Month of Data: January, February, March, April, May, June, July, August, September, October, November, December
# Day of Data: 01, 05, 09, 13, 17, 21, 25, 29
# Time of Data: 00:00, 06:00, 12:00, 18:00
# Pressure Levels: 20 hPa, 30 hPa, 50 hPa, 70 hPa
# Area of Data North: 70°, West: -20°, South: -30°, East: 40° 
# Data Format: NetCDF

PRESSURE_LEVEL = 50.0 # hPa Available levels: 20 hPa, 30 hPa, 50 hPa, 70 hPa
FILEPATH = "/Users/kaimurtagh/Downloads/b00885fac425ddae6f2421ac6486807a.nc"
LEVEL_TOLERANCE = 0.5 # hPa
PERCENTILE = 0.99  # Example: 99th percentile


def load_dataset(filepath):
    try:
        ds = xr.open_dataset(filepath)
    except FileNotFoundError:
        print(f"Error: File not found at {filepath}")
        return None
    except Exception as e:
        print(f"Error loading dataset: {e}")
        return None
    ds['wind_speed'] = np.sqrt(ds['u']**2 + ds['v']**2)
    return ds

def get_wind_summary_stats_for_plotting(ds, lat, lon, requested_level, level_tolerance):
    percentiles_to_calc = [0.95] # ADJUST PERCENTILE HERE!
    speed_p95 = np.nan
    avg_direction = np.nan
    mean_u_val, mean_v_val = np.nan, np.nan

    try:
        subset_speed_da = ds['wind_speed'].sel(
            latitude=lat, longitude=lon, pressure_level=requested_level, method='nearest'
        )
        actual_plev_speed = subset_speed_da.pressure_level.item()

        u_components_da = ds['u'].sel(
            latitude=lat, longitude=lon, pressure_level=requested_level, method='nearest'
        )
        actual_plev_u = u_components_da.pressure_level.item()

        v_components_da = ds['v'].sel(
            latitude=lat, longitude=lon, pressure_level=requested_level, method='nearest'
        )
        actual_plev_v = v_components_da.pressure_level.item()

    except KeyError:
        # print(f"Debug: KeyError during sel for {lat=}, {lon=}, {requested_level=}")
        return {"speed_p95": np.nan, "average_wind_direction": np.nan, "mean_u": np.nan, "mean_v": np.nan}
    except Exception as e: # Catch other potential sel errors
        # print(f"Debug: Exception during sel for {lat=}, {lon=}, {requested_level=}: {e}")
        return {"speed_p95": np.nan, "average_wind_direction": np.nan, "mean_u": np.nan, "mean_v": np.nan}


    mismatch_plev = (abs(actual_plev_speed - requested_level) > level_tolerance or
                     abs(actual_plev_u - requested_level) > level_tolerance or
                     abs(actual_plev_v - requested_level) > level_tolerance)

    if mismatch_plev:
        # print(f"Warning at ({lat:.2f},{lon:.2f}): Plev mismatch. Req: {requested_level}, "
        #       f"Actual_spd: {actual_plev_speed:.2f}. Data ignored.") # Simplified warning
        pass # Data remains NaN
    else:
        if subset_speed_da.size > 0 and not subset_speed_da.isnull().all():
            speed_p95 = float(np.nanpercentile(subset_speed_da.values, percentiles_to_calc[0] * 100))

        if u_components_da.size > 0 and not u_components_da.isnull().all() and \
           v_components_da.size > 0 and not v_components_da.isnull().all():
            
            mean_u_val = u_components_da.mean(skipna=True).item()
            mean_v_val = v_components_da.mean(skipna=True).item()

            if not (np.isnan(mean_u_val) or np.isnan(mean_v_val)):
                if not (mean_u_val == 0 and mean_v_val == 0):
                    avg_direction = (np.arctan2(mean_u_val, mean_v_val) * 180 / np.pi + 360) % 360
    
    return {"speed_p95": speed_p95, "average_wind_direction": avg_direction, "mean_u": mean_u_val, "mean_v": mean_v_val}


def plot_wind_map_eu(ds, requested_level, eu_bounds, grid_resolution, level_tolerance, percentile_to_calc):
    lats_plot_grid = np.arange(eu_bounds['min_lat'], eu_bounds['max_lat'] + grid_resolution[0], grid_resolution[0])
    lons_plot_grid = np.arange(eu_bounds['min_lon'], eu_bounds['max_lon'] + grid_resolution[1], grid_resolution[1])

    speed_data_p95 = np.full((len(lats_plot_grid), len(lons_plot_grid)), np.nan)
    u_component_data = np.full((len(lats_plot_grid), len(lons_plot_grid)), np.nan)
    v_component_data = np.full((len(lats_plot_grid), len(lons_plot_grid)), np.nan)

    # Global warning for pressure level (as before)
    available_plevs = ds.pressure_level.values
    min_diff_to_available = np.min(np.abs(available_plevs - requested_level))
    if min_diff_to_available > level_tolerance * 5:
        print(f"Global Warning: Requested pressure level {requested_level} hPa is significantly different "
              f"from any available levels. Map results may be sparse or all NaN.")

    print(f"Calculating wind statistics for {len(lats_plot_grid) * len(lons_plot_grid)} grid points at {requested_level} hPa...")
    for i, lat_val in enumerate(tqdm(lats_plot_grid, desc="Latitudes")):
        for j, lon_val in enumerate(lons_plot_grid):
            stats = get_wind_summary_stats_for_plotting(ds, lat_val, lon_val, requested_level, level_tolerance)
            
            speed_data_p95[i, j] = stats['speed_p95']
            u_component_data[i, j] = stats['mean_u']
            v_component_data[i, j] = stats['mean_v']

    # --- NEW: Clip data to actual dataset boundaries ---
    actual_data_min_lat = ds.latitude.min().item()
    actual_data_max_lat = ds.latitude.max().item()
    actual_data_min_lon = ds.longitude.min().item()
    actual_data_max_lon = ds.longitude.max().item()
    
    print(f"Dataset geographical boundaries: Lat ({actual_data_min_lat:.2f} to {actual_data_max_lat:.2f}), "
          f"Lon ({actual_data_min_lon:.2f} to {actual_data_max_lon:.2f})")

    # Create 2D grids of the plot coordinates
    lon_mesh_plot, lat_mesh_plot = np.meshgrid(lons_plot_grid, lats_plot_grid)

    # Create a mask for points in the plot grid that are OUTSIDE the data's actual extent
    buffer = grid_resolution[0] / 2  # Example buffer, can be adjusted
    outside_data_mask = (
        (lat_mesh_plot < actual_data_min_lat - buffer) | (lat_mesh_plot > actual_data_max_lat + buffer) |
        (lon_mesh_plot < actual_data_min_lon - buffer) | (lon_mesh_plot > actual_data_max_lon + buffer)
    )

    speed_data_p95[outside_data_mask] = np.nan
    u_component_data[outside_data_mask] = np.nan
    v_component_data[outside_data_mask] = np.nan

    print("Generating map...")
    fig = plt.figure(figsize=(12, 10))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([eu_bounds['min_lon'] - 2, eu_bounds['max_lon'] + 2, 
                   eu_bounds['min_lat'] - 2, eu_bounds['max_lat'] + 2], crs=ccrs.PlateCarree())

    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='aliceblue')
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    
    # Colormap: make NaNs transparent
    current_cmap = plt.cm.get_cmap('viridis').copy()
    current_cmap.set_bad(color='none')

    valid_speeds = speed_data_p95[~np.isnan(speed_data_p95)]
    vmin, vmax = (np.min(valid_speeds), np.max(valid_speeds)) if len(valid_speeds) > 0 else (0, 1)

    if vmin == vmax and len(valid_speeds) > 0: vmin -= 0.5; vmax += 0.5
    elif vmin == vmax and len(valid_speeds) == 0: vmin = 0; vmax = 1
    
    mesh = ax.pcolormesh(lon_mesh_plot, lat_mesh_plot, speed_data_p95, 
                         transform=ccrs.PlateCarree(), 
                         cmap=current_cmap, 
                         vmin=vmin, vmax=vmax, shading='auto')
    plt.colorbar(mesh, ax=ax, orientation='vertical', label=f'{percentile_to_calc*100}th Percentile Wind Speed (m/s) at {requested_level} hPa')

    K = 2 
    if not (np.all(np.isnan(u_component_data)) or np.all(np.isnan(v_component_data))):
        valid_u = u_component_data[~np.isnan(u_component_data)]
        valid_v = v_component_data[~np.isnan(v_component_data)]
        
        max_abs_u = np.max(np.abs(valid_u)) if len(valid_u) > 0 else 1.0
        max_abs_v = np.max(np.abs(valid_v)) if len(valid_v) > 0 else 1.0
        quiver_scale_val = max(max_abs_u, max_abs_v, 1.0) * 20 

        ax.quiver(lon_mesh_plot[::K, ::K], lat_mesh_plot[::K, ::K], 
                  u_component_data[::K, ::K], v_component_data[::K, ::K],
                  transform=ccrs.PlateCarree(), 
                  color='red', scale=quiver_scale_val, width=0.003, headwidth=3, headlength=5)
    else:
        print("Warning: No valid U/V component data to plot wind direction arrows (possibly all outside data bounds or NaN).")

    ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
    plt.title(f'{percentile_to_calc*100}th Percentile Wind Speed and Average Direction - EU - {requested_level} hPa (Data clipped to source extent)')
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Generate a wind map for the EU.")
    parser.add_argument(
        "filepath", type=str, nargs='?', default=FILEPATH, help="Path to the NetCDF file"
    )
    parser.add_argument("--level", type=float, default=PRESSURE_LEVEL, help="Pressure level in hPa")
    parser.add_argument("--level_tolerance", type=float, default=LEVEL_TOLERANCE, help="Tolerance for pressure level matching (hPa)")
    parser.add_argument("--min_lat", type=float, default=30.0, help="Minimum latitude for plot extent")
    parser.add_argument("--max_lat", type=float, default=80.0, help="Maximum latitude for plot extent")
    parser.add_argument("--min_lon", type=float, default=-20.0, help="Minimum longitude for plot extent")
    parser.add_argument("--max_lon", type=float, default=40.0, help="Maximum longitude for plot extent")
    parser.add_argument("--lat_step", type=float, default=1, help="Latitude step for calculation grid")
    parser.add_argument("--lon_step", type=float, default=1, help="Longitude step for calculation grid")
    # Remove the default value for percentile
    parser.add_argument("--percentile", type=float, help="Percentile to calculate (e.g., 0.99 for 99th percentile)")
    args = parser.parse_args()

    ds = load_dataset(args.filepath)
    if ds is None:
        return

    if 'latitude' not in ds.coords or 'longitude' not in ds.coords:
        print("Error: Dataset does not contain 'latitude' or 'longitude' coordinates.")
        return
    if 'pressure_level' not in ds.coords:
        print("Error: Dataset does not contain 'pressure_level' coordinate.")
        return

    eu_bounds_config = {
        'min_lat': args.min_lat, 'max_lat': args.max_lat,
        'min_lon': args.min_lon, 'max_lon': args.max_lon
    }
    grid_res_config = (args.lat_step, args.lon_step)

    # Use the global PERCENTILE variable if args.percentile is not provided
    percentile_to_use = args.percentile if args.percentile is not None else PERCENTILE

    # Pass the percentile argument to the function
    plot_wind_map_eu(ds, args.level, eu_bounds_config, grid_res_config, args.level_tolerance, percentile_to_use)

if __name__ == "__main__":
    main()