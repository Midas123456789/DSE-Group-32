from LH2_conventional import Conventional_LH
from LH2_flying_wing import FlyingWing_LH
from LH2_tandem import Tandem_LH
import numpy as np
import matplotlib.pyplot as plt

# Prepare storage
altitudes = np.linspace(12000, 18000, 3)

results = {
    "Flying Wing": {"b": [], "V": [], "P_required": [], "wing_area": [], "M_total": []},
    "Tandem": {"b": [], "V": [], "P_required": [], "wing_area": [], "M_total": []},
    "Conventional": {"b": [], "V": [], "P_required": [], "wing_area": [], "M_total": []}
}

# Run models at each altitude
for alt in altitudes:
    for name, cls in [("Flying Wing", FlyingWing_LH), ("Tandem", Tandem_LH), ("Conventional", Conventional_LH)]:
        model = cls(altitude=alt)
        sol = model.solve()
        results[name]["b"].append(sol["b"])
        results[name]["V"].append(sol["V"])
        results[name]["P_required"].append(sol["P_required"])
        results[name]["wing_area"].append(sol["wing_area"])
        results[name]["M_total"].append(sol["M_total"])

# Plotting
variables = {
    "b": "Wingspan [m]",
    "V": "Cruise Speed [m/s]",
    "P_required": "Power Required [W]",
    "wing_area": "Wing Area [m²]",
    "M_total": "Total Mass [kg]"
}

fig, axs = plt.subplots(len(variables), 1, figsize=(8, 15), sharex=True)

for i, (key, ylabel) in enumerate(variables.items()):
    ax = axs[i]
    for config, data in results.items():
        ax.plot(altitudes, data[key], label=config)
    ax.set_ylabel(ylabel)
    ax.grid(True)
    if i == 0:
        ax.set_title("Performance Metrics vs Altitude")
    if i == len(variables) - 1:
        ax.set_xlabel("Altitude [m]")
    ax.legend()

plt.tight_layout()
plt.show()
