import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import interp1d

# --- User-defined inputs ---
FILEPATH = "/Users/kaimurtagh/Downloads/8e9011e42c73e534ff7c9200a183f98.nc" # Replace with your file path

PROCESSING_MODE = "grid"  # Options: "single" or "grid"

# Inputs for "single" mode
LATITUDE_SINGLE = 53
LONGITUDE_SINGLE = -8

# Inputs for "grid" mode
# Smaller ranges are for faster testing.
GRID_BOUNDARIES = {
    "north": 80,  # Example: 80 (Your target)
    "south": 30,  # Example: 30 (Your target)
    "west": -20,   # Example: -20 (Your target)
    "east": 40    # Example: 40 (Your target)
}
GRID_RESOLUTION_LAT = 1  # Degrees Example: 2.5 or 1
GRID_RESOLUTION_LON = 1  # Degrees Example: 2.5 or 1

PERCENTILES_STR = "1,5,10,25,50,75,90,95,99"

# --- COESA 1976 Standard Atmosphere Data ---
COESA_DATA = {
    "altitudes_m": [0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000,
                    15000, 20000, 25000, 30000, 40000, 50000, 60000, 70000, 80000],
    "pressures_hpa": [1013.25, 898.76, 795.00, 701.20, 616.60, 540.48, 472.17, 411.06,
                      356.51, 308.00, 265.00, 120.40, 54.75, 25.10, 11.97, 2.63, 0.55,
                      0.11, 0.02, 0.00]
}

# Create interpolation function once globally
COESA_INTERP_FUNC = interp1d(COESA_DATA["pressures_hpa"], COESA_DATA["altitudes_m"],
                             bounds_error=False, fill_value=np.nan)

def pressure_to_altitude_coesa(pressure_hpa):
    """Converts pressure (hPa) to altitude (m) using the 1976 COESA model."""
    if pressure_hpa <= 0 or np.isnan(pressure_hpa):
        return np.nan
    return COESA_INTERP_FUNC(pressure_hpa)

# --- Data Loading and Processing Functions ---
def load_dataset(filepath):
    """Loads the NetCDF dataset and calculates wind speed."""
    try:
        ds = xr.open_dataset(filepath)
    except FileNotFoundError:
        print(f"Error: File not found at {filepath}")
        return None
    except Exception as e:
        print(f"Error loading dataset: {e}")
        return None

    if 'u' not in ds or 'v' not in ds:
        print("Error: Dataset must contain 'u' and 'v' wind components.")
        return None

    ds['wind_speed'] = np.sqrt(ds['u']**2 + ds['v']**2)
    ds['wind_speed'].attrs['units'] = 'm/s'
    ds['wind_speed'].attrs['long_name'] = 'Wind Speed'

    required_coords = ['latitude', 'longitude', 'pressure_level']
    for coord in required_coords:
        if coord not in ds.coords:
            if coord in ds.data_vars: # If it's a data variable, promote it
                ds = ds.set_coords(coord)
            else:
                print(f"Error: Required coordinate '{coord}' not found in the dataset. Available coords: {list(ds.coords.keys())}, Available vars: {list(ds.data_vars.keys())}")
                return None
    return ds

def get_wind_profile_at_point(ds_point_data, percentiles_to_calc, pressure_levels_for_profile):
    """
    Calculates wind speed percentiles for a single selected data point.
    ds_point_data: xarray.Dataset or DataArray already subsetted for a single lat/lon.
    pressure_levels_for_profile: 1D array of pressure levels to use for the profile.
    """
    altitudes_m = np.array([pressure_to_altitude_coesa(p) for p in pressure_levels_for_profile])
    valid_alt_indices = ~np.isnan(altitudes_m)
    altitudes_m_valid = altitudes_m[valid_alt_indices]
    pressure_levels_for_processing_valid = pressure_levels_for_profile[valid_alt_indices]

    sort_indices = np.argsort(altitudes_m_valid)
    altitudes_m_sorted = altitudes_m_valid[sort_indices]
    final_pressure_levels_to_iterate = pressure_levels_for_processing_valid[sort_indices]

    wind_speed_percentiles_at_alt = {p: [] for p in percentiles_to_calc}
    
    point_pressure_levels = ds_point_data['pressure_level'].values if 'pressure_level' in ds_point_data.coords else []

    for p_level_hpa in final_pressure_levels_to_iterate:
        if p_level_hpa not in point_pressure_levels:
            for percentile_val in percentiles_to_calc:
                wind_speed_percentiles_at_alt[percentile_val].append(np.nan)
            continue
        
        try:
            wind_speed_series_at_plev = ds_point_data['wind_speed'].sel(pressure_level=p_level_hpa).values
            wind_speed_values = wind_speed_series_at_plev.flatten() # Handle multiple time points, if any
            wind_speed_values = wind_speed_values[~np.isnan(wind_speed_values)]

            if wind_speed_values.size == 0:
                for percentile_val in percentiles_to_calc:
                    wind_speed_percentiles_at_alt[percentile_val].append(np.nan)
            else:
                for percentile_val in percentiles_to_calc:
                    calculated_percentile = np.percentile(wind_speed_values, percentile_val)
                    wind_speed_percentiles_at_alt[percentile_val].append(calculated_percentile)
        except Exception: # Catches errors if p_level_hpa not found or other issues
            # print(f"Warning: Error processing p_level {p_level_hpa} for point. Appending NaNs. Error: {e}")
            for percentile_val in percentiles_to_calc:
                wind_speed_percentiles_at_alt[percentile_val].append(np.nan)

    for p_val in percentiles_to_calc:
        wind_speed_percentiles_at_alt[p_val] = np.array(wind_speed_percentiles_at_alt[p_val])
    
    return altitudes_m_sorted, wind_speed_percentiles_at_alt

def get_wind_profile_percentiles_single(ds, target_lat, target_lon, percentiles_to_calc):
    """Calculates wind speed percentiles for a single location."""
    try:
        data_at_loc = ds.sel(latitude=target_lat, longitude=target_lon, method='nearest')
    except Exception as e:
        print(f"Error selecting location ({target_lat}, {target_lon}): {e}")
        return None, None, None

    # Use pressure levels from the selected data for this specific point
    if 'pressure_level' not in data_at_loc.coords:
        print(f"Error: 'pressure_level' not found as a coordinate in the selected data for {target_lat}, {target_lon}.")
        return None, None, data_at_loc # Return data_at_loc for title context
        
    pressure_levels_hpa = data_at_loc['pressure_level'].values
    
    altitudes_m_sorted, wind_speed_percentiles_at_alt = get_wind_profile_at_point(
        data_at_loc, percentiles_to_calc, pressure_levels_hpa
    )
    return altitudes_m_sorted, wind_speed_percentiles_at_alt, data_at_loc

def get_average_wind_profile_percentiles_grid(ds, grid_bounds, lat_res, lon_res, percentiles_to_calc):
    """Calculates average wind speed percentiles over a grid."""
    lats_req = np.arange(grid_bounds['south'], grid_bounds['north'] + lat_res, lat_res)
    lons_req = np.arange(grid_bounds['west'], grid_bounds['east'] + lon_res, lon_res)

    dataset_lats = ds.latitude.values
    dataset_lons = ds.longitude.values
    
    # Filter requested lats/lons to be within dataset's actual coordinate range
    lats = lats_req[(lats_req >= dataset_lats.min()) & (lats_req <= dataset_lats.max())]
    lons = lons_req[(lons_req >= dataset_lons.min()) & (lons_req <= dataset_lons.max())]

    if len(lats) == 0 or len(lons) == 0:
        print("Error: No valid grid points found within the dataset's coordinate range based on requested grid.")
        print(f"Dataset Lat range: {dataset_lats.min():.2f}-{dataset_lats.max():.2f}, Lon range: {dataset_lons.min():.2f}-{dataset_lons.max():.2f}")
        print(f"Requested Lat range for grid: {lats_req.min():.2f}-{lats_req.max():.2f}, Lon range for grid: {lons_req.min():.2f}-{lons_req.max():.2f}")
        return None, None

    grid_points = []
    for lat_g in lats:
        for lon_g in lons:
            grid_points.append((lat_g, lon_g))
    
    if not grid_points:
        print("No grid points to process. Check boundaries and resolution against dataset coverage.")
        return None, None

    print(f"Attempting to process {len(grid_points)} grid points...")

    common_pressure_levels_hpa = ds['pressure_level'].values
    if common_pressure_levels_hpa.ndim > 1: # Should be 1D
        print("Warning: 'pressure_level' coordinate is not 1D. Attempting to use first slice.")
        # This might require dataset-specific handling if pressure levels vary spatially/temporally
        first_slice_indices = {dim: 0 for dim in ds['pressure_level'].dims if dim != 'pressure_level'} # type: ignore
        common_pressure_levels_hpa = ds['pressure_level'].isel(first_slice_indices).values


    all_points_percentile_data = {p: [] for p in percentiles_to_calc}
    processed_points_count = 0
    ref_altitudes = None # To store the altitude array from the first successfully processed point

    for i, (lat_g, lon_g) in enumerate(grid_points):
        # print(f"  Processing point {i+1}/{len(grid_points)}: Lat {lat_g:.2f}, Lon {lon_g:.2f}")
        try:
            data_at_loc = ds.sel(latitude=lat_g, longitude=lon_g, method='nearest')
        except Exception as e:
            # print(f"  Skipping point ({lat_g:.2f}, {lon_g:.2f}) due to selection error: {e}")
            continue
        
        altitudes_m_sorted_point, percentile_values_at_alt = get_wind_profile_at_point(
            data_at_loc, percentiles_to_calc, common_pressure_levels_hpa
        )

        if altitudes_m_sorted_point is None or len(altitudes_m_sorted_point) == 0:
            # print(f"  Skipping point ({lat_g:.2f}, {lon_g:.2f}) due to profile calculation error or no valid altitudes.")
            continue
        
        if ref_altitudes is None: # First successful point defines the reference altitude structure
            ref_altitudes = altitudes_m_sorted_point
        # Ensure current point's altitude structure matches the reference
        elif not np.array_equal(ref_altitudes, altitudes_m_sorted_point):
            # This can happen if common_pressure_levels_hpa leads to different valid altitude sets for different points
            # (e.g., some pressure levels entirely missing for a point -> shorter altitude array)
            # For simplicity, we skip points that don't conform.
            # A more complex solution would be to interpolate all profiles to a fixed altitude grid.
            # print(f"  Warning: Altitude array for point ({lat_g:.2f}, {lon_g:.2f}) has a different structure. Skipping.")
            # print(f"    Ref altitudes len: {len(ref_altitudes)}, Point altitudes len: {len(altitudes_m_sorted_point)}")
            continue

        # Validate data lengths before appending
        valid_point_data = True
        current_point_data_for_stacking = {}
        for p_val in percentiles_to_calc:
            if len(percentile_values_at_alt[p_val]) != len(ref_altitudes):
                # print(f"  Warning: Data length mismatch for percentile {p_val} at point ({lat_g:.2f}, {lon_g:.2f}). Expected {len(ref_altitudes)}, got {len(percentile_values_at_alt[p_val])}. Skipping point.")
                valid_point_data = False
                break
            current_point_data_for_stacking[p_val] = percentile_values_at_alt[p_val]
        
        if valid_point_data:
            for p_val in percentiles_to_calc:
                 all_points_percentile_data[p_val].append(current_point_data_for_stacking[p_val])
            processed_points_count += 1
    
    print(f"Successfully processed data for {processed_points_count} out of {len(grid_points)} potential grid points.")
    if processed_points_count == 0:
        print("No data successfully processed from any grid points. Cannot compute average.")
        return None, None

    # Average the collected percentile data
    averaged_percentiles = {p: np.full(len(ref_altitudes) if ref_altitudes is not None else 0, np.nan) for p in percentiles_to_calc} # type: ignore
    if ref_altitudes is None: # Should not happen if processed_points_count > 0
        print("Error: Reference altitudes not set, cannot average.")
        return None, None
        
    for p_val in percentiles_to_calc:
        if not all_points_percentile_data[p_val]: # If list is empty for this percentile
            continue # Will keep the pre-filled NaNs
        
        # Each item in all_points_percentile_data[p_val] is a 1D array of (n_altitudes)
        # Stacking them creates a 2D array: (n_processed_points, n_altitudes)
        stacked_data = np.array(all_points_percentile_data[p_val])
        
        if stacked_data.ndim == 2 and stacked_data.shape[0] > 0 :
             averaged_percentiles[p_val] = np.nanmean(stacked_data, axis=0) # Average along points axis
        elif stacked_data.ndim == 1 and processed_points_count == 1: # Single valid point in grid mode
             averaged_percentiles[p_val] = stacked_data 
        # Else, if stacked_data is problematic, it keeps the NaNs

    return ref_altitudes, averaged_percentiles

def plot_wind_profile(altitudes, percentile_data, percentiles_requested, title_suffix):
    """Plots the wind speed percentile profiles."""
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(10, 8)) # Use ax for plotting

    colors = cm.turbo(np.linspace(0.05, 0.95, len(percentiles_requested)))
    sorted_percentiles_for_plotting = sorted(percentiles_requested)

    if altitudes is None or len(altitudes) == 0:
        print("Error: Cannot plot, altitudes array is missing or empty.")
        ax.text(0.5, 0.5, "No altitude data available for plotting.", 
                horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        plt.title(f'Wind Speed Percentiles {title_suffix} - NO ALTITUDE DATA', fontsize=14)
        plt.show()
        return

    plotted_anything = False
    for i, p_val in enumerate(sorted_percentiles_for_plotting):
        if p_val not in percentile_data or percentile_data[p_val] is None:
            # print(f"No data for percentile {p_val} in plot.")
            continue
        speeds = percentile_data[p_val]
        
        if not isinstance(altitudes, np.ndarray) or altitudes.ndim != 1: continue # Should be 1D
        if not isinstance(speeds, np.ndarray) or speeds.ndim != 1 or speeds.size == 0: continue
        if len(altitudes) != len(speeds):
            # print(f"Plotting warning: Mismatch in lengths of altitudes ({len(altitudes)}) and speeds ({len(speeds)}) for percentile {p_val}%.")
            continue

        valid_indices = ~np.isnan(speeds) & ~np.isnan(altitudes)
        current_altitudes = altitudes[valid_indices]
        current_speeds = speeds[valid_indices]
        
        # Need at least k+1 points for cubic spline (k=3), so at least 4 points.
        if len(current_altitudes) > 3 and len(current_speeds) > 3:
            try:
                interp_func = interp1d(current_altitudes, current_speeds, kind='cubic', bounds_error=False, fill_value="extrapolate")
                smooth_altitudes = np.linspace(current_altitudes.min(), current_altitudes.max(), 300)
                smooth_speeds = interp_func(smooth_altitudes)
                ax.plot(smooth_speeds, smooth_altitudes, 
                         label=f'{p_val}%', color=colors[i], linewidth=2)
                plotted_anything = True
            except ValueError as e: # If interpolation fails (e.g., not enough unique points)
                # print(f"Could not interpolate for percentile {p_val}%: {e}. Plotting raw data points.")
                if len(current_altitudes) > 0: # Still try to plot raw if interpolation failed
                    ax.plot(current_speeds, current_altitudes, 
                            label=f'{p_val}% (raw)', color=colors[i], linewidth=1.5, linestyle='--')
                    plotted_anything = True
        elif len(current_altitudes) > 0 : # If not enough points for cubic, plot raw
             ax.plot(current_speeds, current_altitudes, 
                     label=f'{p_val}% (raw)', color=colors[i], linewidth=1.5, linestyle='--')
             plotted_anything = True

    if not plotted_anything:
        ax.text(0.5, 0.5, "No valid data to plot after filtering.", 
                horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        plt.title(f'Wind Speed Percentiles {title_suffix} - NO VALID DATA TO PLOT', fontsize=14)
    else:
        plt.title(f'Wind Speed Percentiles {title_suffix}', fontsize=14)

    plt.xlabel('Wind Speed [m/s]', fontsize=12)
    plt.ylabel('Altitude [m]', fontsize=12)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    if plotted_anything: # Only show legend if something was plotted
        ax.legend(title="Percentiles", fontsize=10, loc='best')
    plt.tight_layout()
    plt.show()

def main():
    ds = load_dataset(FILEPATH)
    if ds is None:
        return

    try:
        percentiles_to_calculate = [float(p.strip()) for p in PERCENTILES_STR.split(',')]
        if not all(0 <= p <= 100 for p in percentiles_to_calculate): # Allow 0 and 100
            raise ValueError("Percentiles must be between 0 and 100 (inclusive).")
    except ValueError as e:
        print(f"Error parsing percentiles: {e}")
        return

    if PROCESSING_MODE == "single":
        print(f"Processing single point: Lat {LATITUDE_SINGLE}, Lon {LONGITUDE_SINGLE}")
        altitudes, percentile_data, loc_data = get_wind_profile_percentiles_single(
            ds, LATITUDE_SINGLE, LONGITUDE_SINGLE, percentiles_to_calculate
        )
        if altitudes is None or percentile_data is None:
            print("Failed to generate percentile data for single point. Exiting.")
            # Try to plot with whatever title info we have, even if data is bad
            title_lat = LATITUDE_SINGLE
            title_lon = LONGITUDE_SINGLE
            if loc_data is not None: # loc_data might be returned even if altitudes/percentile_data are None
                title_lat = loc_data.latitude.item() if 'latitude' in loc_data.coords else LATITUDE_SINGLE
                title_lon = loc_data.longitude.item() if 'longitude' in loc_data.coords else LONGITUDE_SINGLE
            title_suffix = f'(Lat: {title_lat:.2f}, Lon: {title_lon:.2f})'
            plot_wind_profile(altitudes, percentile_data, percentiles_to_calculate, title_suffix) # Will show "NO DATA"
            return
        
        # Extract actual lat/lon from selected data for title, if available
        title_lat = loc_data.latitude.item() if loc_data is not None and 'latitude' in loc_data.coords else LATITUDE_SINGLE
        title_lon = loc_data.longitude.item() if loc_data is not None and 'longitude' in loc_data.coords else LONGITUDE_SINGLE
        title_suffix = f'(Lat: {title_lat:.2f}, Lon: {title_lon:.2f})'
        plot_wind_profile(altitudes, percentile_data, percentiles_to_calculate, title_suffix)

    elif PROCESSING_MODE == "grid":
        print(f"Processing grid: N:{GRID_BOUNDARIES['north']}, S:{GRID_BOUNDARIES['south']}, W:{GRID_BOUNDARIES['west']}, E:{GRID_BOUNDARIES['east']}")
        print(f"Resolution: Lat {GRID_RESOLUTION_LAT} deg, Lon {GRID_RESOLUTION_LON} deg")
        
        altitudes, avg_percentile_data = get_average_wind_profile_percentiles_grid(
            ds, GRID_BOUNDARIES, GRID_RESOLUTION_LAT, GRID_RESOLUTION_LON, percentiles_to_calculate
        )

        title_suffix = (f'(Grid Avg: {GRID_BOUNDARIES["south"]}°N-{GRID_BOUNDARIES["north"]}°N, '
                        f'{GRID_BOUNDARIES["west"]}°E-{GRID_BOUNDARIES["east"]}°E; '
                        f'{GRID_RESOLUTION_LAT}°x{GRID_RESOLUTION_LON}° res)')
        
        if altitudes is None or avg_percentile_data is None:
            print("Failed to generate averaged percentile data for grid. Exiting.")
            plot_wind_profile(altitudes, avg_percentile_data, percentiles_to_calculate, title_suffix) # Will show "NO DATA"
            return
        
        plot_wind_profile(altitudes, avg_percentile_data, percentiles_to_calculate, title_suffix)

    else:
        print(f"Error: Unknown PROCESSING_MODE '{PROCESSING_MODE}'. Choose 'single' or 'grid'.")

if __name__ == "__main__":
    main()