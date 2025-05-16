#CREATES A WIND MAP FOR THE EU USING HISTORICAL DATA FROM ERA5
# File: b00885fac425ddae6f2421ac6486807a.nc
# Product Type: Reanalysis
# Variables: U-component of wind, V-component of wind
# Year of Data: 2022, 2023, 2024
# Month of Data: January, February, March, April, May, June, July, August, September, October, November, December
# Day of Data: 01, 05, 09, 13, 17, 21, 25, 29
# Time of Data: 00:00, 06:00, 12:00, 18:00
# Pressure Levels: 20 hPa, 30 hPa, 50 hPa, 70 hPa (Example, actual levels from file are used)
# Area of Data North: 70°, West: -20°, South: -30°, East: 40°
# Data Format: NetCDF

import numpy as np
import xarray as xr
import argparse
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from tqdm import tqdm
import matplotlib.colors
import matplotlib # For colormaps

# (Keep your existing file header/comments)

PRESSURE_LEVEL = 50.0 # hPa (Default initial pressure level for the interactive plot)
FILEPATH = "/Users/kaimurtagh/Downloads/b00885fac425ddae6f2421ac6486807a.nc" # Replace with your actual file path
PERCENTILE = 0.99
PRECOMPUTE_PLEV_STEP = 10.0 # hPa step for precomputation

def load_dataset(filepath):
    try:
        ds = xr.open_dataset(filepath)
    except FileNotFoundError:
        print(f"Error: File not found at {filepath}")
        return None
    except Exception as e:
        print(f"Error loading dataset: {e}")
        return None
    if 'pressure_level' in ds.coords:
        ds = ds.sortby('pressure_level')
    ds['wind_speed'] = np.sqrt(ds['u']**2 + ds['v']**2)
    return ds

def get_wind_summary_stats_for_plotting(ds_point_data, requested_level, target_percentile_value):
    """
    Calculates wind statistics for a single lat/lon point by interpolating to the requested_level.
    ds_point_data is expected to be a DataArray/Dataset already selected for a single lat/lon.
    """
    calculated_speed_percentile = np.nan
    avg_direction = np.nan
    mean_u_val, mean_v_val = np.nan, np.nan
    try:
        u_interp = ds_point_data['u'].interp(pressure_level=requested_level)
        v_interp = ds_point_data['v'].interp(pressure_level=requested_level)
        wind_speed_interp = np.sqrt(u_interp**2 + v_interp**2)

        if wind_speed_interp.notnull().any():
            calculated_speed_percentile = float(np.nanpercentile(wind_speed_interp.data, target_percentile_value * 100))
            if u_interp.notnull().any():
                mean_u_val = u_interp.mean(skipna=True).item()
            if v_interp.notnull().any():
                mean_v_val = v_interp.mean(skipna=True).item()
            if not (np.isnan(mean_u_val) or np.isnan(mean_v_val)):
                if not (mean_u_val == 0 and mean_v_val == 0):
                    avg_direction = (np.arctan2(mean_u_val, mean_v_val) * 180 / np.pi + 360) % 360
    except Exception:
        pass
    return {"speed_percentile": calculated_speed_percentile, "average_wind_direction": avg_direction, "mean_u": mean_u_val, "mean_v": mean_v_val}

def precompute_wind_data_for_levels(ds, eu_bounds, grid_resolution, target_percentile_value, plev_step_hz):
    lats_plot_grid = np.arange(eu_bounds['min_lat'], eu_bounds['max_lat'] + grid_resolution[0], grid_resolution[0])
    lons_plot_grid = np.arange(eu_bounds['min_lon'], eu_bounds['max_lon'] + grid_resolution[1], grid_resolution[1])

    available_plevs_ds = ds.pressure_level.values
    ds_plev_min_val, ds_plev_max_val = available_plevs_ds.min(), available_plevs_ds.max()
    
    # Ensure correct order for arange
    start_plev_for_arange = min(ds_plev_min_val, ds_plev_max_val)
    end_plev_for_arange = max(ds_plev_min_val, ds_plev_max_val)

    # Define the pressure levels for precomputation
    precomputation_levels = np.arange(start_plev_for_arange, end_plev_for_arange + plev_step_hz, plev_step_hz)
    # Ensure the exact dataset boundaries are included if they don't fall on a step
    precomputation_levels = np.unique(np.concatenate((precomputation_levels, [ds_plev_min_val, ds_plev_max_val])))
    precomputation_levels.sort() # Ensure sorted

    print(f"\nPrecomputing data for {len(precomputation_levels)} pressure levels:")
    print(f"  Range: {precomputation_levels.min():.1f} to {precomputation_levels.max():.1f} hPa.")
    print(f"  Levels: {precomputation_levels}")

    all_precomputed_data = {}

    # Static geographical mask (calculated once)
    actual_data_min_lat = ds.latitude.min().item()
    actual_data_max_lat = ds.latitude.max().item()
    actual_data_min_lon = ds.longitude.min().item()
    actual_data_max_lon = ds.longitude.max().item()
    lon_mesh_plot_temp, lat_mesh_plot_temp = np.meshgrid(lons_plot_grid, lats_plot_grid)
    buffer_temp = grid_resolution[0] / 2
    outside_data_mask_static = (
        (lat_mesh_plot_temp < actual_data_min_lat - buffer_temp) | (lat_mesh_plot_temp > actual_data_max_lat + buffer_temp) |
        (lon_mesh_plot_temp < actual_data_min_lon - buffer_temp) | (lon_mesh_plot_temp > actual_data_max_lon + buffer_temp)
    )

    for p_level in tqdm(precomputation_levels, desc="Precomputing Levels"):
        speed_grid = np.full((len(lats_plot_grid), len(lons_plot_grid)), np.nan)
        u_grid = np.full((len(lats_plot_grid), len(lons_plot_grid)), np.nan)
        v_grid = np.full((len(lats_plot_grid), len(lons_plot_grid)), np.nan)

        for i, lat_val in enumerate(lats_plot_grid):
            for j, lon_val in enumerate(lons_plot_grid):
                # Select data for this lat/lon point once
                ds_point_data = ds.sel(latitude=lat_val, longitude=lon_val, method='nearest')
                stats = get_wind_summary_stats_for_plotting(ds_point_data, p_level, target_percentile_value)
                speed_grid[i, j] = stats['speed_percentile']
                u_grid[i, j] = stats['mean_u']
                v_grid[i, j] = stats['mean_v']
        
        speed_grid[outside_data_mask_static] = np.nan
        u_grid[outside_data_mask_static] = np.nan
        v_grid[outside_data_mask_static] = np.nan

        all_precomputed_data[p_level] = {
            'speed_grid': speed_grid,
            'u_grid': u_grid,
            'v_grid': v_grid
        }
    print("Precomputation complete.")
    return all_precomputed_data, precomputation_levels


def plot_wind_map_eu_interactive(ds_metadata, precomputed_data_cache, precomputed_plevs_list, 
                                 initial_pressure_level_arg, eu_bounds, grid_resolution, target_percentile_value):
    fig = plt.figure(figsize=(12, 11))
    ax_map = fig.add_axes([0.05, 0.15, 0.9, 0.8], projection=ccrs.PlateCarree())
    ax_slider = fig.add_axes([0.2, 0.05, 0.6, 0.03])

    lats_plot_grid = np.arange(eu_bounds['min_lat'], eu_bounds['max_lat'] + grid_resolution[0], grid_resolution[0])
    lons_plot_grid = np.arange(eu_bounds['min_lon'], eu_bounds['max_lon'] + grid_resolution[1], grid_resolution[1])
    lon_mesh_plot, lat_mesh_plot = np.meshgrid(lons_plot_grid, lats_plot_grid)

    K = 2

    print(f"Dataset geographical boundaries (for plot masking): Lat ({ds_metadata.latitude.min().item():.2f} to {ds_metadata.latitude.max().item():.2f}), "
          f"Lon ({ds_metadata.longitude.min().item():.2f} to {ds_metadata.longitude.max().item():.2f})")

    ax_map.set_extent([eu_bounds['min_lon'] - 2, eu_bounds['max_lon'] + 2,
                       eu_bounds['min_lat'] - 2, eu_bounds['max_lat'] + 2], crs=ccrs.PlateCarree())
    ax_map.add_feature(cfeature.LAND, facecolor='lightgray')
    ax_map.add_feature(cfeature.OCEAN, facecolor='aliceblue')
    ax_map.add_feature(cfeature.COASTLINE)
    ax_map.add_feature(cfeature.BORDERS, linestyle=':')
    gl = ax_map.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
    gl.top_labels = False; gl.right_labels = False

    current_cmap = matplotlib.colormaps.get_cmap('viridis').copy()
    current_cmap.set_bad(color='none')

    initial_speed_data_for_mesh = np.full((len(lats_plot_grid), len(lons_plot_grid)), np.nan)
    mesh = ax_map.pcolormesh(lon_mesh_plot, lat_mesh_plot, initial_speed_data_for_mesh,
                             transform=ccrs.PlateCarree(), cmap=current_cmap, shading='auto')
    cb = fig.colorbar(mesh, ax=ax_map, orientation='vertical', pad=0.02)
    Q = ax_map.quiver(lon_mesh_plot[::K, ::K], lat_mesh_plot[::K, ::K],
                      np.zeros_like(lon_mesh_plot[::K, ::K]), np.zeros_like(lon_mesh_plot[::K, ::K]),
                      transform=ccrs.PlateCarree(), color='red', scale=100, width=0.003, headwidth=3, headlength=5,
                      pivot='middle')

    def update_plot_from_precomputed(slider_plev_value):
        # Find the closest precomputed pressure level to the slider's current value
        # This is important if the slider can output values not exactly in precomputed_plevs_list
        actual_plev_to_use = precomputed_plevs_list[np.argmin(np.abs(precomputed_plevs_list - slider_plev_value))]

        print(f"\nUpdating plot using precomputed data for pressure level: {actual_plev_to_use:.2f} hPa (slider val: {slider_plev_value:.2f})")

        if actual_plev_to_use not in precomputed_data_cache:
            print(f"  Error: Pressure level {actual_plev_to_use:.2f} not found in precomputed_data_cache.")
            return

        data_for_level = precomputed_data_cache[actual_plev_to_use]
        speed_percentile_data_grid = data_for_level['speed_grid']
        u_component_data_grid = data_for_level['u_grid']
        v_component_data_grid = data_for_level['v_grid']

        mesh.set_array(speed_percentile_data_grid)
        valid_speeds = speed_percentile_data_grid[~np.isnan(speed_percentile_data_grid)]
        vmin, vmax = (np.min(valid_speeds), np.max(valid_speeds)) if len(valid_speeds) > 0 else (0, 1)
        if vmin == vmax and len(valid_speeds) > 0: vmin -= 0.5; vmax += 0.5
        elif vmin == vmax and len(valid_speeds) == 0: vmin = 0; vmax = 1
        mesh.set_clim(vmin, vmax)

        cb.set_label(f'{target_percentile_value*100:.0f}th Pctl. Wind Speed (m/s)\nPrecomputed at {actual_plev_to_use:.1f} hPa')
        cb.update_normal(mesh)

        if not (np.all(np.isnan(u_component_data_grid)) or np.all(np.isnan(v_component_data_grid))):
            valid_u = u_component_data_grid[~np.isnan(u_component_data_grid)]
            valid_v = v_component_data_grid[~np.isnan(v_component_data_grid)]
            if len(valid_u) > 0 and len(valid_v) > 0 :
                max_abs_u = np.max(np.abs(valid_u)); max_abs_v = np.max(np.abs(valid_v))
                current_quiver_scale_value = max(max_abs_u, max_abs_v, 0.1) * 20
                if current_quiver_scale_value == 0: current_quiver_scale_value = 100
                Q.set_UVC(u_component_data_grid[::K, ::K], v_component_data_grid[::K, ::K])
                Q.scale = current_quiver_scale_value
            else:
                 Q.set_UVC(np.zeros_like(u_component_data_grid[::K, ::K]), np.zeros_like(v_component_data_grid[::K, ::K])); Q.scale = 100
        else:
            Q.set_UVC(np.zeros_like(u_component_data_grid[::K, ::K]), np.zeros_like(v_component_data_grid[::K, ::K])); Q.scale = 100
            print("  Warning: No valid U/V data in precomputed set for this level.")

        ax_map.set_title(f'{target_percentile_value*100:.0f}th Pctl. Wind Speed & Avg. Direction - EU\nPrecomputed at {actual_plev_to_use:.1f} hPa', fontsize=12)
        fig.canvas.draw_idle()
        print(f"Plot updated for {actual_plev_to_use:.2f} hPa using precomputed data.")

    # Slider setup using the actual precomputed pressure levels
    slider_plev_min, slider_plev_max = precomputed_plevs_list.min(), precomputed_plevs_list.max()
    
    # Find the precomputed level closest to the initial argument for slider init
    initial_slider_val = precomputed_plevs_list[np.argmin(np.abs(precomputed_plevs_list - initial_pressure_level_arg))]

    # Determine a sensible valstep for the slider
    # If precomputed_plevs_list are uniformly PRECOMPUTE_PLEV_STEP apart, this is easy.
    # Otherwise, a small step with "find closest" is robust.
    slider_step = PRECOMPUTE_PLEV_STEP
    if len(precomputed_plevs_list) > 1:
        diffs = np.diff(precomputed_plevs_list)
        # Check if most diffs are close to PRECOMPUTE_PLEV_STEP
        if np.all(np.isclose(diffs, PRECOMPUTE_PLEV_STEP)) or len(diffs) == 0:
            slider_step = PRECOMPUTE_PLEV_STEP
        else: # If not uniform, use a smaller step and rely on "find closest"
            slider_step = max(1.0, np.min(diffs[diffs > 1e-5]) if len(diffs[diffs>1e-5]) > 0 else 1.0)
            print(f"  Note: Precomputed pressure levels are not uniformly {PRECOMPUTE_PLEV_STEP}hPa apart. Slider step adjusted to ~{slider_step:.2f}hPa.")


    pressure_slider = Slider(
        ax=ax_slider,
        label='Pressure (hPa)',
        valmin=slider_plev_min,
        valmax=slider_plev_max,
        valinit=initial_slider_val,
        valstep=slider_step # Slider will step according to this.
                           # update_plot_from_precomputed will find the NEAREST actual precomputed level.
    )
    pressure_slider.on_changed(update_plot_from_precomputed)
    update_plot_from_precomputed(initial_slider_val) # Initial plot draw
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Generate an interactive wind map for the EU with precomputation.")
    parser.add_argument("filepath", type=str, nargs='?', default=FILEPATH, help="Path to the NetCDF file")
    parser.add_argument("--initial_level", type=float, default=PRESSURE_LEVEL, help="Initial pressure level in hPa for the plot")
    parser.add_argument("--min_lat", type=float, default=30.0, help="Minimum latitude")
    parser.add_argument("--max_lat", type=float, default=70.0, help="Maximum latitude")
    parser.add_argument("--min_lon", type=float, default=-25.0, help="Minimum longitude")
    parser.add_argument("--max_lon", type=float, default=40.0, help="Maximum longitude")
    parser.add_argument("--lat_step", type=float, default=1.0, help="Latitude step for grid") # Smaller = slower precomputation
    parser.add_argument("--lon_step", type=float, default=1.0, help="Longitude step for grid") # Smaller = slower precomputation
    parser.add_argument("--percentile", type=float, help=f"Percentile (0-1). Default {PERCENTILE}")
    parser.add_argument("--precompute_step", type=float, default=PRECOMPUTE_PLEV_STEP, help=f"Step for pressure level precomputation (hPa). Default {PRECOMPUTE_PLEV_STEP} hPa.")

    args = parser.parse_args()


    print("Loading dataset...")
    ds = load_dataset(args.filepath)
    if ds is None: return
    for coord in ['latitude', 'longitude', 'pressure_level']:
        if coord not in ds.coords or ds[coord].size == 0:
            print(f"Error: Dataset missing or empty required coordinate: '{coord}'.")
            return

    eu_bounds_config = {'min_lat': args.min_lat, 'max_lat': args.max_lat, 'min_lon': args.min_lon, 'max_lon': args.max_lon}
    grid_res_config = (args.lat_step, args.lon_step)
    percentile_to_use = args.percentile if args.percentile is not None else PERCENTILE
    if not (0 < percentile_to_use < 1):
        print(f"Error: Percentile must be between 0 and 1 (exclusive). Received: {percentile_to_use}"); return

    # Precompute all necessary data
    precomputed_data_cache, precomputed_plevs_list = precompute_wind_data_for_levels(
        ds, eu_bounds_config, grid_res_config, percentile_to_use, args.precompute_step
    )

    if not precomputed_data_cache:
        print("Error: Precomputation resulted in no data. Exiting.")
        return

    # Pass ds for metadata (like actual lat/lon ranges for plot titles/info if needed, and grid setup)
    # but the core data for plotting comes from precomputed_data_cache
    plot_wind_map_eu_interactive(
        ds, precomputed_data_cache, precomputed_plevs_list,
        args.initial_level, eu_bounds_config, grid_res_config, percentile_to_use
    )

if __name__ == "__main__":
    main()

