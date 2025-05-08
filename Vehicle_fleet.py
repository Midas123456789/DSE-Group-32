import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shapely.geometry import Point, box
from math import sqrt


def load_mainland_eu(filepath: str, crs="EPSG:3035"):
    eu = gpd.read_file(filepath)
    eu_members = ['AT', 'BE', 'BG', 'HR', 'CY', 'CZ', 'DK', 'EE', 'FI', 'FR',
                  'DE', 'EL', 'HU', 'IE', 'IT', 'LV', 'LT', 'LU', 'MT', 'NL',
                  'PL', 'PT', 'RO', 'SK', 'SI', 'ES', 'SE']
    eu = eu[eu['CNTR_CODE'].isin(eu_members)].to_crs(crs)
    eu_parts = eu.explode(index_parts=False)
    centroids = eu_parts.geometry.centroid
    mainland = eu_parts[
        (centroids.y.between(900000, 5500000)) &
        (centroids.x.between(2500000, 7500000))
    ]
    return mainland

def create_grid_within_eu(gdf, n_cells=10, crs="EPSG:3035"):
    xmin, ymin, xmax, ymax = gdf.total_bounds
    cell_size = (xmax - xmin) / n_cells
    grid_cells = [
        box(x0, y0, x0 + cell_size, y0 + cell_size)
        for x0 in np.arange(xmin, xmax + cell_size, cell_size)
        for y0 in np.arange(ymin, ymax + cell_size, cell_size)
    ]
    grid = gpd.GeoDataFrame(geometry=grid_cells, crs=crs)
    grid = gpd.sjoin(grid, gdf, how='inner', predicate='intersects')
    grid = grid.drop(columns=[col for col in ['index_right'] if col in grid.columns])
    return grid

def generate_hexagonal_centers(bounds, radius):
    xmin, ymin, xmax, ymax = bounds
    dx = 1.5 * radius
    dy = sqrt(3) * radius
    points = []
    x = xmin
    while x < xmax + radius:
        y_offset = 0 if (int((x - xmin) / dx) % 2 == 0) else dy / 2
        y = ymin + y_offset
        while y < ymax + radius:
            points.append(Point(x, y))
            y += dy
        x += dx
    return points

def create_circles(points, radius, crs="EPSG:3035"):
    gdf = gpd.GeoDataFrame(geometry=points, crs=crs)
    circles = gdf.copy()
    circles["geometry"] = gdf.buffer(radius)
    return gdf, circles

def filter_circles_by_grid(circles, grid):
    circles = gpd.sjoin(circles, grid, how="inner", predicate="intersects")
    circles = circles.drop_duplicates('geometry')
    return circles

def remove_redundant_circles(circles, grid):
    circle_grid_map = gpd.sjoin(circles[['geometry']], grid[['geometry']], how="left", predicate="intersects")
    circle_grid_map['circle_id'] = circle_grid_map.index
    grid_assignments = circle_grid_map.groupby('circle_id').apply(lambda df: set(df.index_right)).to_dict()

    to_remove = set()
    circle_ids = list(grid_assignments.keys())
    for i in range(len(circle_ids)):
        for j in range(len(circle_ids)):
            if i == j or circle_ids[i] in to_remove or circle_ids[j] in to_remove:
                continue
            a, b = circle_ids[i], circle_ids[j]
            set_a, set_b = grid_assignments[a], grid_assignments[b]
            if set_a <= set_b:
                to_remove.add(a if len(set_a) <= len(set_b) else b)
            elif set_b <= set_a:
                to_remove.add(b if len(set_b) <= len(set_a) else a)
    return circles.drop(index=to_remove)

def transform_to_wgs84(*gdfs):
    return [gdf.to_crs(epsg=4326) for gdf in gdfs]

def plot_circles(grid, circles, centers):
    fig, ax = plt.subplots(figsize=(12, 12))
    grid.plot(ax=ax, edgecolor='black', facecolor='none')
    circles.plot(ax=ax, edgecolor='blue', facecolor='none', alpha=0.3, label="150 km Circles")
    centers.plot(ax=ax, color='red', markersize=5, label="Circle Centers")
    plt.title("150 km Circles Covering Mainland EU Grid (Redundant Removed)")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_aspect('equal')
    plt.legend()
    plt.show()


def main(shapefile_path, n_cells=200, radius_m=150_000):
    mainland = load_mainland_eu(shapefile_path)
    grid = create_grid_within_eu(mainland, n_cells=n_cells)
    circle_centers = generate_hexagonal_centers(grid.total_bounds, radius_m)
    centers_gdf, circles_gdf = create_circles(circle_centers, radius_m)
    circles_gdf = filter_circles_by_grid(circles_gdf, grid)
    circles_gdf = remove_redundant_circles(circles_gdf, grid)
    centers_gdf = circles_gdf.copy()
    centers_gdf["geometry"] = centers_gdf["geometry"].centroid
    centers_gdf, circles_gdf, grid = transform_to_wgs84(centers_gdf, circles_gdf, grid)

    plot_circles(grid, circles_gdf, centers_gdf)

    return centers_gdf

def export_circle_coordinates(centers_gdf, output_path=None):
    """
    Converts circle center geometries to a DataFrame with lat/lon and optionally saves to CSV.

    Parameters:
        centers_gdf (GeoDataFrame): Circle centers in EPSG:4326.
        output_path (str, optional): If provided, saves the DataFrame to this CSV path.

    Returns:
        pd.DataFrame: DataFrame with circle center coordinates.
    """
    coords = centers_gdf.geometry.apply(lambda pt: (pt.y, pt.x))  # (lat, lon)
    df_locations = pd.DataFrame(coords.tolist(), columns=["latitude", "longitude"])
    df_locations.index.name = "circle_id"

    if output_path:
        df_locations.to_csv(output_path)

    return df_locations

# Example usage:
centers_gdf = main("Europe_Communes/COMM_RG_01M_2016_3035.shp")
export_circle_coordinates(centers_gdf, "circle_centers.csv")
