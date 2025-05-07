import geopandas as gpd
import matplotlib.pyplot as plt

# Load shapefile
eu = gpd.read_file("Europe_map/CNTR_BN_01M_2024_3035.shp")

# Optional: simplify by selecting only EU countries
# If 'CNTR_CODE' or 'NAME_ENGL' is in the dataframe, you can filter here
# For now, just plot everything

# Plot
# Keep only countries where EU_FLAG is '1'
eu_euonly = eu[eu['EU_FLAG'] == 'T']

# Plot the filtered EU countries
eu_euonly.plot(figsize=(10, 10), color='lightgreen', edgecolor='black')

plt.title("Europe (Eurostat shapefile)")
plt.show()
print(eu.columns)
