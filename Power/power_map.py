from Power.solar_power_model import *  # Must define Power class with average_power_12h()
import plotly.graph_objects as go
import numpy as np

# Define grid
days = np.arange(0, 365)                  # Day of year
latitudes = np.arange(34.5, 81.8)             # Latitude range (rows)

# Compute average power
power_map = np.zeros((len(latitudes), len(days)))  # shape (lat, day)

for i, lat in enumerate(latitudes):
    for j, day in enumerate(days):
        power_model = Power(latitude=lat, day_of_year=day, area=1)
        min_power = 0
        average_power_12h = power_model.average_power_12h()
        #average_power_12h = average_power_12h if average_power_12h > min_power else None
        power_map[i, j] = average_power_12h

# Build figure
fig = go.Figure()

# Add heatmap layer
fig.add_trace(go.Heatmap(
    z=power_map,
    x=days,
    y=latitudes,
    colorscale='Cividis',
    colorbar=dict(title="Avg Power [W]"),
    zsmooth='best',
    showscale=True,
))

# Add contour lines
fig.add_trace(go.Contour(
    z=power_map,
    x=days,
    y=latitudes,
    contours=dict(
        coloring='none',        # Only draw lines
        showlabels=True,        # Show values on lines
        labelfont=dict(size=12, color='black')
    ),
    line=dict(color='black', width=1),
    showscale=False,            # Don't show duplicate colorbar
))

# Layout
fig.update_layout(
    title="Average Solar Power vs Latitude and Day of Year",
    xaxis_title="Day of Year",
    yaxis_title="Latitude [°]",
    template='plotly_white',
)

fig.show()
