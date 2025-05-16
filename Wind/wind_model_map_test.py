# tests/test_wind_analysis.py
import pytest
import numpy as np
import xarray as xr
import pandas as pd
from unittest.mock import patch, MagicMock, ANY

# Import functions and constants from your script
# Assuming wind_analysis.py is in the parent directory, conftest.py handles path.
import wind_model_map
from wind_model_map import (
    load_dataset,
    get_wind_summary_stats_for_plotting,
    plot_wind_map_eu,
    main,
    LEVEL_TOLERANCE as SCRIPT_LEVEL_TOLERANCE, # Alias to avoid name clashes
    PERCENTILE as SCRIPT_PERCENTILE
)

# Helper function to create a mock xarray.Dataset for tests
def create_mock_dataset(
    lat_values=np.array([-10.0, 0.0, 10.0]),
    lon_values=np.array([-5.0, 0.0, 5.0]),
    plev_values=np.array([30.0, 50.0, 70.0]),
    time_values_count=2,
    u_fill_value=10.0,
    v_fill_value=-5.0
):
    """Creates a mock xarray.Dataset with u and v components."""
    times = pd.to_datetime([f'2023-01-01T{h:02d}:00:00' for h in range(time_values_count)])
    
    # Ensure u_fill_value and v_fill_value are broadcastable or match shape
    if isinstance(u_fill_value, (int, float)):
        u_data = np.full((len(times), len(plev_values), len(lat_values), len(lon_values)), float(u_fill_value))
    else: # Assume it's already an array of the correct shape
        u_data = u_fill_value

    if isinstance(v_fill_value, (int, float)):
        v_data = np.full((len(times), len(plev_values), len(lat_values), len(lon_values)), float(v_fill_value))
    else: # Assume it's already an array of the correct shape
        v_data = v_fill_value

    ds = xr.Dataset(
        {
            "u": (("time", "pressure_level", "latitude", "longitude"), u_data),
            "v": (("time", "pressure_level", "latitude", "longitude"), v_data),
        },
        coords={
            "time": times,
            "pressure_level": plev_values,
            "latitude": lat_values,
            "longitude": lon_values,
        },
    )
    return ds

@pytest.fixture
def mock_ds_basic():
    """Pytest fixture for a basic mock dataset."""
    return create_mock_dataset()

@pytest.fixture
def mock_ds_all_nans():
    """Pytest fixture for a mock dataset where u and v are all NaNs."""
    ds = create_mock_dataset(u_fill_value=np.nan, v_fill_value=np.nan)
    # Manually add wind_speed as load_dataset would, but with NaNs
    ds['wind_speed'] = np.sqrt(ds['u']**2 + ds['v']**2) # will be all nans
    return ds


# --- Tests for load_dataset ---
@patch('wind_analysis.xr.open_dataset')
def test_load_dataset_success(mock_open_dataset, mock_ds_basic):
    """Test successful loading and wind_speed calculation."""
    # Prepare a mock dataset without 'wind_speed' initially
    ds_to_load = mock_ds_basic.copy()
    if 'wind_speed' in ds_to_load:
        ds_to_load = ds_to_load.drop_vars('wind_speed')
    mock_open_dataset.return_value = ds_to_load

    loaded_ds = load_dataset("dummy_path.nc")

    mock_open_dataset.assert_called_once_with("dummy_path.nc")
    assert loaded_ds is not None
    assert "wind_speed" in loaded_ds
    expected_wind_speed = np.sqrt(mock_ds_basic['u']**2 + mock_ds_basic['v']**2)
    xr.testing.assert_allclose(loaded_ds['wind_speed'], expected_wind_speed)

@patch('wind_analysis.xr.open_dataset')
@patch('builtins.print')
def test_load_dataset_file_not_found(mock_print, mock_open_dataset):
    """Test load_dataset when file is not found."""
    mock_open_dataset.side_effect = FileNotFoundError("File not found")
    
    loaded_ds = load_dataset("non_existent_path.nc")
    
    assert loaded_ds is None
    mock_print.assert_called_with("Error: File not found at non_existent_path.nc")

@patch('wind_analysis.xr.open_dataset')
@patch('builtins.print')
def test_load_dataset_other_error(mock_print, mock_open_dataset):
    """Test load_dataset with a generic error during opening."""
    mock_open_dataset.side_effect = Exception("Some other error")

    loaded_ds = load_dataset("dummy_path.nc")

    assert loaded_ds is None
    mock_print.assert_called_with("Error loading dataset: Some other error")


# --- Tests for get_wind_summary_stats_for_plotting ---
def test_get_wind_summary_stats_valid_data(mock_ds_basic):
    """Test with valid data and matching pressure level."""
    ds = mock_ds_basic.copy()
    ds['wind_speed'] = np.sqrt(ds['u']**2 + ds['v']**2)
    
    lat, lon = 0.0, 0.0
    requested_level = 50.0
    level_tolerance = 0.5

    # Expected values based on create_mock_dataset defaults (u=10, v=-5)
    # and the hardcoded percentile [0.95] in the function
    expected_speed = np.sqrt(10**2 + (-5)**2) # approx 11.18
    # arctan2(u,v) for meteorological wind direction (from North, clockwise)
    # arctan2(10, -5) gives angle in radians. Convert to degrees.
    # Then (angle * 180/pi + 360) % 360. Note: u,v in arctan2 for wind direction can vary.
    # The code uses arctan2(mean_u, mean_v). If u is east, v is north.
    # arctan2(u,v) typically gives angle from positive x-axis (East) counter-clockwise.
    # (np.arctan2(10, -5) * 180 / np.pi + 360) % 360 = (116.56 + 360) % 360 = 116.56 deg
    # This seems to be mathematical angle (East=0, North=90).
    # Wind direction is usually "from where the wind blows".
    # A wind with u=10 (eastward), v=-5 (southward) comes from NW.
    # The formula (np.arctan2(mean_u_val, mean_v_val) * 180 / np.pi + 360) % 360
    # If u=positive (East), v=positive (North) -> angle in [0, 90] (from SW)
    # If u=10, v=-5. mean_u=10, mean_v=-5. atan2(10, -5) = 2.034 rad = 116.56 deg.
    # This is angle CCW from positive V axis if atan2(U,V). If atan2(Y,X) -> CCW from positive X.
    # Python's atan2(y,x). So it's atan2(mean_u_val, mean_v_val).
    # If u is y-like (meridional) and v is x-like (zonal), then it is standard angle from V.
    # Assuming standard (u=zonal, v=meridional), then atan2(v,u) is common for angle.
    # Script uses atan2(mean_u, mean_v). Let's assume u=x, v=y for this.
    # So, atan2(10, -5) means y=10, x=-5. Angle from negative x-axis.
    # Let's stick to what the code does: (np.arctan2(10, -5) * 180 / np.pi + 360) % 360
    # np.arctan2(10, -5) approx 2.0344 radians. 2.0344 * 180 / np.pi = 116.565 degrees.
    expected_direction = (np.arctan2(10, -5) * 180 / np.pi + 360) % 360

    stats = get_wind_summary_stats_for_plotting(ds, lat, lon, requested_level, level_tolerance)

    assert np.isclose(stats['speed_p95'], expected_speed) # P95 of constant is constant
    assert np.isclose(stats['average_wind_direction'], expected_direction)
    assert np.isclose(stats['mean_u'], 10.0)
    assert np.isclose(stats['mean_v'], -5.0)

def test_get_wind_summary_stats_plev_mismatch(mock_ds_basic):
    """Test when requested pressure level is outside tolerance."""
    ds = mock_ds_basic.copy()
    ds['wind_speed'] = np.sqrt(ds['u']**2 + ds['v']**2)

    stats = get_wind_summary_stats_for_plotting(ds, 0.0, 0.0, 55.0, SCRIPT_LEVEL_TOLERANCE) # 50 is nearest, 55-50=5 > 0.5

    assert np.isnan(stats['speed_p95'])
    assert np.isnan(stats['average_wind_direction'])
    assert np.isnan(stats['mean_u'])
    assert np.isnan(stats['mean_v'])

def test_get_wind_summary_stats_key_error(mock_ds_basic):
    """Test when a key (e.g., 'pressure_level' in selection) is missing, causing KeyError."""
    # Create a dataset that will cause a KeyError by not having 'pressure_level' on subset
    # This is tricky to set up directly, will mock the .sel() behavior instead
    ds_mock = MagicMock(spec=xr.Dataset)
    # Make ds_mock['wind_speed'].sel raise KeyError
    ds_mock.__getitem__.return_value.sel.side_effect = KeyError("Mocked KeyError")

    stats = get_wind_summary_stats_for_plotting(ds_mock, 0.0, 0.0, 50.0, 0.5)
    
    assert np.isnan(stats['speed_p95'])
    assert np.isnan(stats['average_wind_direction'])
    assert np.isnan(stats['mean_u'])
    assert np.isnan(stats['mean_v'])

def test_get_wind_summary_stats_all_nans_input(mock_ds_all_nans):
    """Test with input data that is all NaNs."""
    stats = get_wind_summary_stats_for_plotting(mock_ds_all_nans, 0.0, 0.0, 50.0, 0.5)

    assert np.isnan(stats['speed_p95'])
    assert np.isnan(stats['average_wind_direction'])
    assert np.isnan(stats['mean_u']) # .mean(skipna=True) of all nans is nan
    assert np.isnan(stats['mean_v'])

def test_get_wind_summary_stats_zero_wind(mock_ds_basic):
    """Test case where mean u and v are zero."""
    ds = create_mock_dataset(u_fill_value=0.0, v_fill_value=0.0)
    ds['wind_speed'] = np.sqrt(ds['u']**2 + ds['v']**2)

    stats = get_wind_summary_stats_for_plotting(ds, 0.0, 0.0, 50.0, 0.5)

    assert np.isclose(stats['speed_p95'], 0.0)
    assert np.isnan(stats['average_wind_direction']) # Direction is undefined for zero wind
    assert np.isclose(stats['mean_u'], 0.0)
    assert np.isclose(stats['mean_v'], 0.0)

# --- Tests for plot_wind_map_eu ---
# These tests will mock plotting functions and check data preparation logic.
@patch('wind_analysis.plt') # Mocks matplotlib.pyplot
@patch('wind_analysis.cfeature') # Mocks cartopy.feature
@patch('wind_analysis.ccrs') # Mocks cartopy.crs
@patch('wind_analysis.get_wind_summary_stats_for_plotting')
@patch('builtins.print') # To capture print statements
def test_plot_wind_map_eu_basic_run(mock_print, mock_get_stats, mock_ccrs, mock_cfeature, mock_plt, mock_ds_basic):
    """Test that plot_wind_map_eu runs, calls helpers, and tries to plot."""
    ds = mock_ds_basic.copy()
    ds.attrs['latitude_min'] = ds.latitude.min().item() # Add attributes plot_wind_map_eu expects
    ds.attrs['latitude_max'] = ds.latitude.max().item()
    ds.attrs['longitude_min'] = ds.longitude.min().item()
    ds.attrs['longitude_max'] = ds.longitude.max().item()

    # Mock return value for get_wind_summary_stats_for_plotting
    mock_get_stats.return_value = {
        "speed_p95": 10.0, 
        "average_wind_direction": 90.0, 
        "mean_u": 0.0, 
        "mean_v": 10.0
    }

    requested_level = 50.0
    eu_bounds = {'min_lat': -5.0, 'max_lat': 5.0, 'min_lon': -2.0, 'max_lon': 2.0}
    grid_resolution = (5.0, 2.0) # lat_step, lon_step
    level_tolerance = 0.5
    percentile_to_calc = 0.99 # This is for labels, actual data is P95

    plot_wind_map_eu(ds, requested_level, eu_bounds, grid_resolution, level_tolerance, percentile_to_calc)

    # Check if get_wind_summary_stats_for_plotting was called
    # Lats: -5, 0, 5 (3 points). Lons: -2, 0, 2 (3 points). Total 3*3 = 9 calls.
    assert mock_get_stats.call_count == 9 
    
    # Check if plotting functions were called
    mock_plt.figure.assert_called_once()
    mock_plt.axes.assert_called_once()
    ax_mock = mock_plt.axes.return_value
    ax_mock.set_extent.assert_called_once()
    ax_mock.add_feature.assert_called()
    ax_mock.pcolormesh.assert_called_once()
    mock_plt.colorbar.assert_called_once()
    ax_mock.quiver.assert_called_once() # Assuming some valid U/V data
    mock_plt.title.assert_called_once()
    mock_plt.show.assert_called_once()

    # Check that the title uses the passed percentile_to_calc
    title_call_args = mock_plt.title.call_args[0][0]
    assert f'{percentile_to_calc*100}th Percentile Wind Speed' in title_call_args
    
    # Check global warning for pressure level (if applicable)
    # For this test, requested_level 50.0 is in mock_ds_basic.plev_values, so no warning.
    # To test warning: change requested_level to something far off, e.g., 500.0
    # Example: any( "Global Warning: Requested pressure level" in call[0][0] for call in mock_print.call_args_list)


@patch('wind_analysis.plt')
@patch('wind_analysis.get_wind_summary_stats_for_plotting')
@patch('builtins.print')
def test_plot_wind_map_eu_data_clipping(mock_print, mock_get_stats, mock_plt, mock_ds_basic):
    """Test the data clipping logic based on dataset's actual extent."""
    # Dataset covers a small area
    ds_small_extent = create_mock_dataset(
        lat_values=np.array([0.0]), lon_values=np.array([0.0]), plev_values=np.array([50.0])
    )
    # plot_wind_map_eu reads these from ds directly, not attrs
    # ds_small_extent.attrs['latitude_min'] = 0.0
    # ds_small_extent.attrs['latitude_max'] = 0.0
    # ds_small_extent.attrs['longitude_min'] = 0.0
    # ds_small_extent.attrs['longitude_max'] = 0.0
    
    # Mock stats to return non-NaN, so any NaNs are from clipping
    mock_get_stats.return_value = {"speed_p95": 5.0, "mean_u": 1.0, "mean_v": 1.0, "average_wind_direction": 45.0}

    # Plot bounds are much larger than data extent
    eu_bounds = {'min_lat': -10.0, 'max_lat': 10.0, 'min_lon': -10.0, 'max_lon': 10.0}
    grid_resolution = (5.0, 5.0) # Grid points at -10, -5, 0, 5, 10 for lat/lon

    plot_wind_map_eu(ds_small_extent, 50.0, eu_bounds, grid_resolution, 0.5, 0.99)

    # Get the data passed to pcolormesh
    ax_mock = mock_plt.axes.return_value
    pcolormesh_call_args = ax_mock.pcolormesh.call_args
    # args are (lon_mesh, lat_mesh, speed_data, ...)
    speed_data_plot = pcolormesh_call_args[0][2] # This is speed_data_p95 after clipping

    # Expected: speed_data_plot should have NaNs for points outside ds_small_extent (0,0)
    # Lats for plot: -10, -5, 0, 5, 10. Lons for plot: -10, -5, 0, 5, 10.
    # Center point (index 2,2) corresponds to lat=0, lon=0, which is in data extent.
    # All other points should be NaN due to clipping.
    # The buffer in clipping logic might affect exact number of NaNs,
    # but most points should be NaN.
    
    # Check a point expected to be NaN and one expected to be data
    # (0,0) is lat_idx=2, lon_idx=2 in the 5x5 grid
    # A point far from data, e.g., (-10, -10) -> (idx 0,0) should be NaN
    assert np.isnan(speed_data_plot[0, 0]), "Point far from data should be NaN after clipping"
    # The point (0,0) should be the value from mock_get_stats (5.0)
    # Assuming grid_resolution / 2 buffer doesn't exclude the exact match point.
    assert np.isclose(speed_data_plot[2, 2], 5.0), "Point within data extent should have value"
    
    # Count NaNs: all but the center (or a few points around center due to buffer) should be NaN
    assert np.sum(np.isnan(speed_data_plot)) >= (speed_data_plot.size - 4), "Most points should be NaN due to clipping"


# --- Tests for main ---
@patch('wind_analysis.argparse.ArgumentParser')
@patch('wind_analysis.load_dataset')
@patch('wind_analysis.plot_wind_map_eu')
@patch('builtins.print')
def test_main_successful_run(mock_print, mock_plot_map, mock_load, mock_argparse, mock_ds_basic):
    """Test main function happy path."""
    # Mock argparse
    mock_args = MagicMock()
    mock_args.filepath = "dummy.nc"
    mock_args.level = 50.0
    mock_args.level_tolerance = 0.5
    mock_args.min_lat = 30.0; mock_args.max_lat = 60.0
    mock_args.min_lon = -10.0; mock_args.max_lon = 20.0
    mock_args.lat_step = 1.0; mock_args.lon_step = 1.0
    mock_args.percentile = 0.99 # Explicitly provide percentile

    mock_parser_instance = mock_argparse.return_value
    mock_parser_instance.parse_args.return_value = mock_args

    # Mock load_dataset
    mock_load.return_value = mock_ds_basic # Simulate successful load

    main()

    mock_load.assert_called_once_with(mock_args.filepath)
    mock_plot_map.assert_called_once_with(
        mock_ds_basic,
        mock_args.level,
        ANY, # eu_bounds_config dict
        (mock_args.lat_step, mock_args.lon_step),
        mock_args.level_tolerance,
        mock_args.percentile # Should use the parsed percentile
    )
    # Check eu_bounds_config structure (part of ANY)
    plot_call_args = mock_plot_map.call_args[0]
    eu_bounds_arg = plot_call_args[2]
    assert eu_bounds_arg['min_lat'] == mock_args.min_lat

@patch('wind_analysis.argparse.ArgumentParser')
@patch('wind_analysis.load_dataset')
@patch('wind_analysis.plot_wind_map_eu')
@patch('builtins.print')
def test_main_load_dataset_fails(mock_print, mock_plot_map, mock_load, mock_argparse):
    """Test main when load_dataset returns None."""
    mock_args = MagicMock()
    mock_args.filepath = "fail.nc"
    # ... other args ...
    mock_parser_instance = mock_argparse.return_value
    mock_parser_instance.parse_args.return_value = mock_args

    mock_load.return_value = None # Simulate load failure

    main()

    mock_load.assert_called_once_with(mock_args.filepath)
    mock_plot_map.assert_not_called() # Plotting should not occur

@patch('wind_analysis.argparse.ArgumentParser')
@patch('wind_analysis.load_dataset')
@patch('wind_analysis.plot_wind_map_eu')
@patch('builtins.print')
def test_main_uses_default_percentile(mock_print, mock_plot_map, mock_load, mock_argparse, mock_ds_basic):
    """Test main uses global PERCENTILE when --percentile arg is not given."""
    mock_args = MagicMock()
    mock_args.filepath = "dummy.nc"
    mock_args.level = 50.0
    mock_args.level_tolerance = 0.5
    mock_args.min_lat = 30.0; mock_args.max_lat = 60.0
    mock_args.min_lon = -10.0; mock_args.max_lon = 20.0
    mock_args.lat_step = 1.0; mock_args.lon_step = 1.0
    mock_args.percentile = None # Simulate --percentile not being provided

    mock_parser_instance = mock_argparse.return_value
    mock_parser_instance.parse_args.return_value = mock_args
    mock_load.return_value = mock_ds_basic

    main()

    mock_plot_map.assert_called_once()
    # Check that the SCRIPT_PERCENTILE (the global one) was used
    assert mock_plot_map.call_args[0][5] == SCRIPT_PERCENTILE


@patch('wind_analysis.argparse.ArgumentParser')
@patch('wind_analysis.load_dataset')
@patch('wind_analysis.plot_wind_map_eu')
@patch('builtins.print')
def test_main_missing_coordinates_error(mock_print, mock_plot_map, mock_load, mock_argparse):
    """Test main function error handling for missing coordinates in dataset."""
    mock_args = MagicMock()
    mock_args.filepath = "dummy.nc"
    # ... other default args ...
    mock_args.percentile = 0.99
    mock_parser_instance = mock_argparse.return_value
    mock_parser_instance.parse_args.return_value = mock_args

    # Create a dataset missing 'latitude'
    ds_missing_coords = xr.Dataset({'u': (('time', 'longitude'), np.random.rand(2,2))},
                                   coords={'time': pd.to_datetime(['2023-01-01']),
                                           'longitude': [0,1]})
    mock_load.return_value = ds_missing_coords

    main()

    mock_load.assert_called_once_with(mock_args.filepath)
    mock_print.assert_any_call("Error: Dataset does not contain 'latitude' or 'longitude' coordinates.")
    mock_plot_map.assert_not_called()