import aerosandbox as asb
import aerosandbox.numpy as np
import matplotlib.pyplot as plt
from aerosandbox import Atmosphere

# Constants
W = 4000  # Weight [N]
Cd = 0.04
S = 100   # Wing area [m^2]
V = 40    # Flight speed [m/s]

# Arrays to store results
thrust_values = np.arange(0.05 * W, 1.0 * W, 100)  # Thrust range
altitudes = np.arange(0, 30001, 1000)  # Altitude range
results = []

for T in thrust_values:
    total_time = 0
    total_energy = 0
    for h in range(0, 18000, 1000):
        h_next = h + 1000
        altitude = (h + h_next) / 2  # Midpoint for atmosphere calc
        atmo = Atmosphere(altitude)
        rho = atmo.density()

        D = 0.5 * rho * V**2 * S * Cd
        RC = (T - D) / W  # Rate of climb [m/s]
        if RC <= 0:
            total_time = np.nan
            total_energy = np.nan
            break
        segment_time = 1000 / RC  # Time for 1000 m climb [s]
        total_time += segment_time
        total_energy += segment_time * T  # Energy = Thrust * time

    if not np.isnan(total_time):
        results.append((T, total_time / 3600, total_energy))  # Time in hours

# Convert to arrays for plotting
results = np.array(results)
T_vals = results[:, 0]
times_hr = results[:, 1]
energies = results[:, 2]

# Plotting
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

ax1.plot(T_vals, times_hr, 'b-', label='Time to climb (hr)')
ax2.plot(T_vals, energies / 1e6, 'r--', label='Energy required (MJ)')  # Convert J to MJ

ax1.set_xlabel('Thrust (N)')
ax1.set_ylabel('Time to climb (hr)', color='b')
ax2.set_ylabel('Energy required (MJ)', color='r')
ax1.grid(True)

plt.title('Climb Performance vs Thrust')
fig.tight_layout()
plt.show()
