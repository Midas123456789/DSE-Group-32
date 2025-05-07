import numpy as np
import xarray as xr
import argparse

# Configurable options (MONTH is no longer used as a default for filtering)
LATITUDE = 53.0
LONGITUDE = 30.0
PRESSURE_LEVEL = 70.0
FILEPATH = "/Users/kaimurtagh/Downloads/b00885fac425ddae6f2421ac6486807a.nc" # Replace with a valid path if needed

def load_dataset(filepath):
    """Load a NetCDF dataset and compute wind speed."""
    ds = xr.open_dataset(filepath)
    ds['wind_speed'] = np.sqrt(ds['u']**2 + ds['v']**2)
    # Instantaneous direction might be useful for other analyses, not directly used for the average here
    ds['wind_direction_instantaneous'] = (np.arctan2(ds['u'], ds['v']) * 180 / np.pi + 360) % 360
    return ds

def get_wind_summary_stats(ds, lat, lon, level, speed_percentiles=[0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999]):
    """
    Query wind speed percentiles and average wind direction at a specific
    lat/lon/level, using all available time data.

    Parameters:
    - ds: xarray.Dataset containing u, v, and wind_speed
    - lat: Latitude (float)
    - lon: Longitude (float)
    - level: Pressure level in hPa (float)
    - speed_percentiles: List of percentiles to compute for wind speed

    Returns:
    - Dictionary of requested percentiles for wind speed and the average wind direction.
    """
    # --- Wind Speed Percentiles (using all available time data) ---
    # Select data for the specific point across all time steps
    subset_speed = ds['wind_speed'].sel(
        latitude=lat,
        longitude=lon,
        pressure_level=level,
        method='nearest'
    )

    if subset_speed.size == 0 or subset_speed.isnull().all(): # Check if empty or all NaN
        print(f"Warning: No valid speed data found for lat={lat}, lon={lon}, level={level} across all times. Speed percentiles will be NaN.")
        speed_results = {
            f"speed_p{int(p * 100)}": np.nan
            for p in speed_percentiles
        }
    else:
        speed_results = {
            f"speed_p{int(p * 100)}": float(np.nanpercentile(subset_speed.values, p * 100))
            for p in speed_percentiles
        }

    # --- Average Wind Direction (using all available time data) ---
    # Select u and v components for the specific point across all time steps
    u_components = ds['u'].sel(
        latitude=lat,
        longitude=lon,
        pressure_level=level,
        method='nearest'
    )
    v_components = ds['v'].sel(
        latitude=lat,
        longitude=lon,
        pressure_level=level,
        method='nearest'
    )

    if u_components.size == 0 or v_components.size == 0 or \
       u_components.isnull().all() or v_components.isnull().all():
        print(f"Warning: No valid u/v data found for lat={lat}, lon={lon}, level={level} across all times. Average direction will be NaN.")
        avg_direction = np.nan
    else:
        # Calculate mean of u and v components over all available time
        # xarray's .mean() automatically skips NaNs for numeric types
        mean_u = u_components.mean(skipna=True).item() # .item() converts 0-dim array to scalar
        mean_v = v_components.mean(skipna=True).item()

        if np.isnan(mean_u) or np.isnan(mean_v): # If all values were NaN after selection
            avg_direction = np.nan
            print(f"Note: All u/v components were NaN for lat={lat}, lon={lon}, level={level} across all times. Average direction is NaN.")
        elif mean_u == 0 and mean_v == 0: # Handle calm conditions
            avg_direction = np.nan # Direction is undefined for calm
            print(f"Note: Average wind components are (0,0) (calm) for lat={lat}, lon={lon}, level={level} across all times. Direction is undefined.")
        else:
            # Calculate average direction from mean components (meteorological convention)
            avg_direction = (np.arctan2(mean_u, mean_v) * 180 / np.pi + 360) % 360
    
    direction_result = {"average_wind_direction": float(avg_direction)}

    return {**speed_results, **direction_result}

def validate_inputs(ds, lat, lon, level):
    """
    Validate if the given latitude, longitude, and pressure level exist in the dataset.

    Parameters:
    - ds: xarray.Dataset
    - lat: Latitude (float)
    - lon: Longitude (float)
    - level: Pressure level in hPa (float)

    Raises:
    - ValueError: If any of the inputs are not within the dataset's range.
    """
    # Check latitude
    if lat < ds.latitude.min().item() or lat > ds.latitude.max().item():
        raise ValueError(f"Error: Latitude {lat} is out of bounds. Dataset range: {ds.latitude.min().item()} to {ds.latitude.max().item()}.")

    # Check longitude
    if lon < ds.longitude.min().item() or lon > ds.longitude.max().item():
        raise ValueError(f"Error: Longitude {lon} is out of bounds. Dataset range: {ds.longitude.min().item()} to {ds.longitude.max().item()}.")

    # Check pressure level
    if level not in ds.pressure_level.values:
        raise ValueError(f"Error: Pressure level {level} hPa is not available in the dataset. Available levels: {ds.pressure_level.values}.")

def main():
    parser = argparse.ArgumentParser(description="Query wind speed percentiles and average direction from NetCDF data, using all available time steps.")
    parser.add_argument(
        "filepath",
        type=str,
        help="Path to the NetCDF file",
        nargs='?',
        default=FILEPATH
    )
    parser.add_argument("--lat", type=float, required=False, help="Latitude", default=LATITUDE)
    parser.add_argument("--lon", type=float, required=False, help="Longitude", default=LONGITUDE)
    parser.add_argument("--level", type=float, required=False, help="Pressure level in hPa", default=PRESSURE_LEVEL)
    args = parser.parse_args()

    try:
        ds = load_dataset(args.filepath)
    except FileNotFoundError:
        print(f"Error: File not found at {args.filepath}")
        print("Please provide a valid path to a NetCDF file.")
        return
    except Exception as e:
        print(f"Error loading dataset: {e}")
        return

    try:
        # Validate the inputs
        validate_inputs(ds, args.lat, args.lon, args.level)
    except ValueError as e:
        print(e)
        return

    # Call without the month argument
    result = get_wind_summary_stats(ds, args.lat, args.lon, args.level)

    print(f"\nWind speed percentiles and average direction at lat={args.lat}, lon={args.lon}, level={args.level} hPa (using all available data):")
    for k, v in result.items():
        if isinstance(v, float):
            print(f"  {k}: {v:.2f}")
        else:
            print(f"  {k}: {v}")

if __name__ == "__main__":
    main()