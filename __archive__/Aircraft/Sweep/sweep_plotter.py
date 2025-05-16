import csv
import matplotlib.pyplot as plt

# === Read data from CSV === #
filename = "solar_sweep_results.csv"
results = []

with open(filename, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        results.append({
            "Irradiance": float(row["Irradiance"]),
            "Span": float(row["Span"]),
            "Wing Area": float(row["Wing Area"]),
            "Battery Mass": float(row["Battery Mass"]),
            "Skin Mass": float(row["Skin Mass"]),
            "Spar Mass": float(row["Spar Mass"]),
            "Structural Mass": float(row["Structural Mass"]),
            "Total Weight": float(row["Total Weight"]),
            "Velocity": float(row["Velocity"]),
            "Wing Loading": float(row["Wing Loading"]),
        })

# === Extract data for plotting === #
irr_vals = [r["Irradiance"] for r in results]
span_vals = [r["Span"] for r in results]
area_vals = [r["Wing Area"] for r in results]
bat_mass_vals = [r["Battery Mass"] for r in results]
skin_mass_vals = [r["Skin Mass"] for r in results]
spar_mass_vals = [r["Spar Mass"] for r in results]
struct_mass_vals = [r["Structural Mass"] for r in results]
weight_vals = [r["Total Weight"] for r in results]
velocity_vals = [r["Velocity"] for r in results]
wing_loading_vals = [r["Wing Loading"] for r in results]

# === Plotting === #
plt.figure(figsize=(12, 18))

def make_subplot(index, x, y, ylabel, title):
    ax = plt.subplot(5, 2, index)
    ax.plot(x, y, marker='o')
    ax.set_xlabel("Solar Irradiance (W/m²)")
    ax.set_ylabel(ylabel)
    ax.grid(True)
    # Place the title inside the plot (top-left corner)
    x_pos = x[0] + 0.02 * (x[-1] - x[0])
    y_pos = max(y) - 0.05 * (max(y) - min(y))
    ax.text(x_pos, y_pos, title, fontsize=10, fontweight='bold',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

make_subplot(1, irr_vals, span_vals, "Wingspan (m)", "Wingspan vs Irradiance")
make_subplot(2, irr_vals, area_vals, "Wing Area (m²)", "Wing Area vs Irradiance")
make_subplot(3, irr_vals, bat_mass_vals, "Battery Mass (kg)", "Battery Mass vs Irradiance")
make_subplot(4, irr_vals, skin_mass_vals, "Skin Mass (kg)", "Skin Mass vs Irradiance")
make_subplot(5, irr_vals, spar_mass_vals, "Spar Mass (kg)", "Spar Mass vs Irradiance")
make_subplot(6, irr_vals, struct_mass_vals, "Structural Mass (kg)", "Structural Mass vs Irradiance")
make_subplot(7, irr_vals, weight_vals, "Total Weight (N)", "Total Weight vs Irradiance")
make_subplot(8, irr_vals, velocity_vals, "Velocity (m/s)", "Velocity vs Irradiance")
make_subplot(9, irr_vals, wing_loading_vals, "Wing Loading (N/m²)", "Wing Loading vs Irradiance")

plt.tight_layout()
plt.show()
