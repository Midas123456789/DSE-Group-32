import geopandas as gpd
from shapely.geometry import Point
import matplotlib.pyplot as plt
import math

# Load the EU shapefile
eu = gpd.read_file("europe_10km.shp")

# Ensure the coordinate reference system is in meters for accurate distance calculations
eu = eu.to_crs(epsg=3035)  # ETRS89 / LAEA Europe

fig, ax = plt.subplots(figsize=(10, 10))
eu.plot(ax=ax, color='lightgray', edgecolor='black')
plt.title("Map of Europe")
plt.tight_layout()
plt.show(block=True)
print("Test")
def generate_hex_grid(polygon, sensor_range):
    # Calculate horizontal and vertical spacing between drone positions
    dx = 1.5 * sensor_range
    dy = math.sqrt(3) * sensor_range

    # Get the bounding box of the polygon
    minx, miny, maxx, maxy = polygon.bounds

    # Generate grid points
    points = []
    y = miny
    row = 0
    while y < maxy:
        x = minx + (0.75 * sensor_range if row % 2 else 0)
        while x < maxx:
            point = Point(x, y)
            if polygon.contains(point):
                points.append(point)
            x += dx
        y += dy
        row += 1
    return points

# Assuming the EU shapefile has a single geometry; if multiple, you may need to merge them
eu_polygon = eu.geometry.union_all()

# Define sensor range in meters
sensor_range = 50000  # 50 km

# Generate drone positions
drone_positions = generate_hex_grid(eu_polygon, sensor_range)

# Create a GeoDataFrame for drone positions
drone_gdf = gpd.GeoDataFrame(geometry=drone_positions, crs=eu.crs)
print(f"Number of drone positions: {len(drone_positions)}")

# Plot the EU boundaries and drone positions
fig, ax = plt.subplots(figsize=(10, 10))
eu.plot(ax=ax, color='white', edgecolor='black')
drone_gdf.plot(ax=ax, color='blue', markersize=5)
plt.title(f"Drone Coverage Over EU with Sensor Range {sensor_range/1000} km")
plt.show()

# def compute_drone_positions(area_size, sensor_range):
#     positions = []
#
#     # Step size for hexagonal tiling
#     dx = 1.5 * sensor_range
#     dy = math.sqrt(3) * sensor_range
#
#     y = 0
#     row = 0
#     while y < area_size:
#         x_offset = 0 if row % 2 == 0 else 0.75 * sensor_range
#         x = x_offset
#         while x < area_size:
#             positions.append((x, y))
#             x += dx
#         y += dy
#         row += 1
#
#     return positions
#
# def plot_coverage(area_size, positions, sensor_range):
#     fig, ax = plt.subplots()
#     ax.set_xlim(0, area_size)
#     ax.set_ylim(0, area_size)
#     ax.set_aspect('equal')
#
#     for x, y in positions:
#         circle = plt.Circle((x, y), sensor_range, color='blue', alpha=0.3)
#         ax.add_patch(circle)
#         ax.plot(x, y, 'ko')  # Drone position
#
#     plt.title(f"Drone Coverage with {len(positions)} Drones")
#     plt.grid(True)
#     plt.show()
#
# # Example use
# area_size = 100  # e.g., 100x100 meters
# sensor_range = 10  # e.g., 10 meter coverage radius
#
# positions = compute_drone_positions(area_size, sensor_range)
# plot_coverage(area_size, positions, sensor_range)
#
# print(f"Total drones needed: {len(positions)}")