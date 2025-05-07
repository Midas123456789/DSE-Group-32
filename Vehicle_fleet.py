import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Point, box
from math import sqrt

def create_grid_within_eu(gdf, n_cells=10, crs="EPSG:3035"):
    xmin, ymin, xmax, ymax = gdf.total_bounds
    cell_size = (xmax - xmin) / n_cells

    grid_cells = []
    for x0 in np.arange(xmin, xmax + cell_size, cell_size):
        for y0 in np.arange(ymin, ymax + cell_size, cell_size):
            x1 = x0 + cell_size
            y1 = y0 + cell_size
            grid_cells.append(box(x0, y0, x1, y1))

    grid = gpd.GeoDataFrame(geometry=grid_cells, crs=crs)
    grid = gpd.sjoin(grid, gdf, how='inner', predicate='intersects')
    grid = grid.drop(columns=[col for col in ['index_right'] if col in grid.columns])
    return grid

# Load EU shapefile and filter
eu = gpd.read_file("Europe_Communes/COMM_RG_01M_2016_3035.shp")
eu_members = ['AT', 'BE', 'BG', 'HR', 'CY', 'CZ', 'DK', 'EE', 'FI', 'FR',
              'DE', 'EL', 'HU', 'IE', 'IT', 'LV', 'LT', 'LU', 'MT', 'NL',
              'PL', 'PT', 'RO', 'SK', 'SI', 'ES', 'SE']
eu_euonly = eu[eu['CNTR_CODE'].isin(eu_members)].to_crs(epsg=3035)
eu_parts = eu_euonly.explode(index_parts=False)

# Filter mainland parts by approximate bounding box
centroids = eu_parts.geometry.centroid
mainland_parts = eu_parts[
    (centroids.y.between(900000, 5500000)) &  # y is northing
    (centroids.x.between(2500000, 7500000))   # x is easting
]

# Create grid
grid = create_grid_within_eu(mainland_parts, n_cells=200, crs="EPSG:3035")

# Radius in meters
radius_m = 150_000  # 150 km

# Hexagonal packing
dx = 1.5 * radius_m
dy = sqrt(3) * radius_m

xmin, ymin, xmax, ymax = grid.total_bounds
circle_points = []
x = xmin
while x < xmax + radius_m:
    y_offset = 0 if (int((x - xmin) / dx) % 2 == 0) else dy / 2
    y = ymin + y_offset
    while y < ymax + radius_m:
        circle_points.append(Point(x, y))
        y += dy
    x += dx

# Create circles
centers_gdf = gpd.GeoDataFrame(geometry=circle_points, crs="EPSG:3035")
circles_gdf = centers_gdf.copy()
circles_gdf["geometry"] = centers_gdf.buffer(radius_m)

# Spatial join to filter only circles that intersect grid
circles_gdf = gpd.sjoin(circles_gdf, grid, how="inner", predicate="intersects")
circles_gdf = circles_gdf.drop_duplicates('geometry')
centers_gdf = circles_gdf.copy()
centers_gdf["geometry"] = centers_gdf["geometry"].centroid

centers_gdf_wgs84 = centers_gdf.to_crs(epsg=4326)

# Print longitude and latitude of each center
for idx, row in centers_gdf_wgs84.iterrows():
    lon, lat = row.geometry.x, row.geometry.y
    print(f"Longitude: {lon:.6f}, Latitude: {lat:.6f}")

# Plot
fig, ax = plt.subplots(figsize=(12, 12))
grid.plot(ax=ax, edgecolor='black', facecolor='none')
circles_gdf.plot(ax=ax, edgecolor='blue', facecolor='none', alpha=0.3, label="150 km Circles")
centers_gdf.plot(ax=ax, color='red', markersize=5, label="Circle Centers")
plt.title("150 km Circles Covering Mainland EU Grid")
ax.set_xlabel("Easting (m)")
ax.set_ylabel("Northing (m)")
ax.set_aspect('equal')
plt.legend()
plt.show()
