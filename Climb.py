import aerosandbox as asb
import aerosandbox.numpy as np
import matplotlib.pyplot as plt
from aerosandbox import Atmosphere

# Constants
W = 10000  # Weight [N]
D = 400   # Drag [N]
V = 50  # [m/s]
battery_capacity_Wh = 200000  # Battery capacity [Wh] total
solar_panel_power_W = D*V*24/8   # Solar panel input [W] including night (average)
print("Solar panel power (W):", solar_panel_power_W)



thrust_values = np.arange( D, 1.0 * W, 50)  
results = []

for T in thrust_values:
    total_time = 0
    total_energy = 0
    feasible = True

    for h in range(0, 18000, 1000):
        h_next = h + 1000
        altitude = (h + h_next) / 2
        atmo = Atmosphere(altitude)
        rho = atmo.density()
        RC = (T * V - D * V) / W  # Rate of climb [m/s]

        if RC <= 0:
            feasible = False
            break

        segment_time = 1000 / RC  # Time for 1000 m climb [s]
        total_time += segment_time
        total_energy += segment_time * (T*V)

    if feasible:
        # Now compute available energy from solar + battery
        total_time_hr = total_time / 3600
        solar_energy_Wh = solar_panel_power_W * total_time_hr
        total_available_energy_Wh = battery_capacity_Wh + solar_energy_Wh
        total_available_energy_MJ = total_available_energy_Wh * 3600 / 1e6

        results.append((T, total_time_hr, total_energy, total_available_energy_MJ))

results = np.array(results)
T_vals = results[:, 0]
times_hr = results[:, 1]
energies = results[:, 2]
available_energies = results[:, 3]


fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

ax1.plot(T_vals, times_hr, 'b-', label='Time to climb (hr)')
ax2.plot(T_vals, energies / 1e6, 'r--', label='Energy required for climb flight(MJ)')
ax2.plot(T_vals, available_energies, 'g:', label='Energy Available (MJ)') 


ax1.set_xlabel('Thrust (N)')
ax1.set_ylabel('Time to climb (hr)', color='b')
ax2.set_ylabel('Energy (MJ)', color='r')
ax1.grid(True)
plt.title('Climb Performance vs Thrust')


lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc='upper right')

fig.tight_layout()
plt.show()
