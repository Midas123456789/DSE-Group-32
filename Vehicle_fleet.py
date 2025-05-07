import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np


def create_grid_within_eu(gdf, n_cells=10, crs="EPSG:4326"):
    """Create square grid that covers the EU region"""
    import shapely

    # Get the bounds of the EU region
    xmin, ymin, xmax, ymax = gdf.total_bounds

    # Calculate the cell size (square grid)
    cell_size = (xmax - xmin) / n_cells

    # Create grid cells
    grid_cells = []
    for x0 in np.arange(xmin, xmax + cell_size, cell_size):
        for y0 in np.arange(ymin, ymax + cell_size, cell_size):
            x1 = x0 + cell_size
            y1 = y0 + cell_size
            poly = shapely.geometry.box(x0, y0, x1, y1)
            grid_cells.append(poly)

    # Create GeoDataFrame for the grid
    grid = gpd.GeoDataFrame(grid_cells, columns=['geometry'], crs=crs)

    # Perform spatial join to keep only the grid cells that intersect with the EU boundary
    grid = grid.sjoin(gdf, how='inner').drop_duplicates('geometry')

    return grid


# Load shapefile
eu = gpd.read_file("Europe_Communes/COMM_RG_01M_2016_3035.shp")
eu_members = ['AT', 'BE', 'BG', 'HR', 'CY', 'CZ', 'DK', 'EE', 'FI', 'FR',
              'DE', 'EL', 'HU', 'IE', 'IT', 'LV', 'LT', 'LU', 'MT', 'NL',
              'PL', 'PT', 'RO', 'SK', 'SI', 'ES', 'SE']

# Filter the GeoDataFrame
eu_euonly = eu[eu['CNTR_CODE'].isin(eu_members)]

# Convert to WGS84 for filtering
eu_euonly = eu_euonly.to_crs(epsg=4326)

# Explode multipolygons into individual polygons
eu_parts = eu_euonly.explode(index_parts=False)

# Filter individual polygons by bounding box (approximate mainland Europe)
mainland_parts = eu_parts[
    eu_parts.geometry.centroid.y.between(34, 72) &
    eu_parts.geometry.centroid.x.between(-15, 35)
    ]

# Create grid for the mainland EU region
grid = create_grid_within_eu(mainland_parts, n_cells=200, crs="EPSG:4326")

# Plot only the grid in green with proper lat/lon scale
fig, ax = plt.subplots(figsize=(10, 10))
grid.plot(ax=ax, edgecolor='black', facecolor='lightgreen', linestyle='--', alpha=0.7)
plt.title("Mainland EU Grid (Green)")
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")
ax.set_aspect('equal')  # Ensures 1:1 scale for lat/lon
plt.show()

