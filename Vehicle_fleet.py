import geopandas as gpd
import matplotlib.pyplot as plt

# Load shapefile
eu = gpd.read_file("Europe_map/CNTR_BN_01M_2024_3035.shp")

# Filter only EU member countries
eu_euonly = eu[eu['EU_FLAG'] == 'T']

# Convert to WGS84 for filtering
eu_euonly = eu_euonly.to_crs(epsg=4326)

# Explode multipolygons into individual polygons
eu_parts = eu_euonly.explode(index_parts=False)

# Filter individual polygons by bounding box (approximate mainland Europe)
mainland_parts = eu_parts[
    eu_parts.geometry.centroid.y.between(34, 72) &
    eu_parts.geometry.centroid.x.between(-15, 35)
]


# Plot result
mainland_parts.plot(figsize=(10, 10), color='lightblue', edgecolor='black')
plt.title("Mainland EU (exploded and filtered)")
plt.show()
