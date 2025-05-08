import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
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


def create_grid(gdf, n_cells=10, crs="EPSG:3035"):
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


def compute_circle_centers(bounds, radius):
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


def create_circle_geometries(points, radius, crs="EPSG:3035"):
    gdf = gpd.GeoDataFrame(geometry=points, crs=crs)
    circles = gdf.copy()
    circles["geometry"] = gdf.buffer(radius)
    return gdf, circles


def filter_circles(circles, grid):
    circles = gpd.sjoin(circles, grid, how="inner", predicate="intersects")
    circles = circles.drop_duplicates('geometry')
    return circles


def remove_redundancies(circles, grid):
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


def to_wgs84(*gdfs):
    return [gdf.to_crs(epsg=4326) for gdf in gdfs]


def plot_all(grid, circles, centers):
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


def generate_coverage_data(shapefile_path, n_cells=200, radius_m=150_000):
    mainland = load_mainland_eu(shapefile_path)
    grid = create_grid(mainland, n_cells=n_cells)
    circle_points = compute_circle_centers(grid.total_bounds, radius_m)
    centers_gdf, circles_gdf = create_circle_geometries(circle_points, radius_m)
    circles_gdf = filter_circles(circles_gdf, grid)
    circles_gdf = remove_redundancies(circles_gdf, grid)

    centers_gdf = circles_gdf.copy()
    centers_gdf["geometry"] = centers_gdf["geometry"].centroid

    centers_gdf, circles_gdf, grid = to_wgs84(centers_gdf, circles_gdf, grid)
    return grid, circles_gdf, centers_gdf

def export_coverage_data_shp(grid, circles_gdf, centers_gdf, folder_name="coverage_shapefiles"):
    """
    Saves shapefiles of grid, circles, and centers into a folder in the main directory.

    Parameters:
        grid (GeoDataFrame): Grid polygons.
        circles_gdf (GeoDataFrame): Coverage circles.
        centers_gdf (GeoDataFrame): Circle centers.
        folder_name (str): Name of the output folder in the main directory.
    """
    main_dir = os.getcwd()
    output_dir = os.path.join(main_dir, folder_name)
    os.makedirs(output_dir, exist_ok=True)

    grid.to_file(os.path.join(output_dir, "grid.shp"))
    circles_gdf.to_file(os.path.join(output_dir, "circles.shp"))
    centers_gdf.to_file(os.path.join(output_dir, "centers.shp"))

    print(f"Exported shapefiles to: {output_dir}")


def plot_coverage(folder_name="coverage_shapefiles"):
    """
    Loads exported shapefiles and plots grid, circles, and centers.

    Parameters:
        folder_name (str): Folder containing grid.shp, circles.shp, and centers.shp.
    """
    main_dir = os.getcwd()
    shp_dir = os.path.join(main_dir, folder_name)

    grid = gpd.read_file(os.path.join(shp_dir, "grid.shp"))
    circles = gpd.read_file(os.path.join(shp_dir, "circles.shp"))
    centers = gpd.read_file(os.path.join(shp_dir, "centers.shp"))

    plot_all(grid, circles, centers)


# Generate and export
grid, circles, centers = generate_coverage_data("Europe_Communes/COMM_RG_01M_2016_3035.shp")
export_coverage_data_shp(grid, circles, centers)

# Later, or in another script: just plot from shapefiles
plot_coverage()