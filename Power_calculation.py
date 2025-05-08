import aerosandbox as asb
import aerosandbox.numpy as np
from aerosandbox import cas
from numpy import pi
from aerosandbox import Atmosphere
import matplotlib.pyplot as plt

# Constants
altitude = 20000  # meters
D = 400           # Thrust [N]
v_o = 40          # Flight speed [m/s]
atmo = Atmosphere(altitude)
rho = atmo.density()

# Sweep upper bound on propeller radius from 0.75 to 4 m
r_max_vals = np.arange(0.75, 6, 0.05)

# Storage
r_opt_list = []
rpm_list = []
pitch_list = []
torque_list = []
power_list = []
rpm_torque_ratio_list = []

# Loop over upper bound values of r
for r_max in r_max_vals:
    opti = asb.Opti()

    # Optimization variables
    r = opti.variable(init_guess=1, lower_bound=0.5, upper_bound=r_max)  # m
    RPM = opti.variable(init_guess=600, lower_bound=500, upper_bound=20000)  # RPM

    # Derived quantities
    Adisk = pi * r**2
    T = D
    power = 0.5 * rho * T * v_o * (cas.sqrt(T / (0.5 * rho * Adisk * v_o**2) + 1) + 1)
    pitch = v_o / (RPM / 60)

    # Objective: minimize power
    opti.minimize(power)

    # Solve
    sol = opti.solve()

    # Extract optimized values
    r_opt = sol.value(r)
    RPM_opt = sol.value(RPM)
    P_opt = sol.value(power)
    pitch_opt = sol.value(pitch)
    omega = 2 * pi * RPM_opt / 60
    torque_opt = P_opt / omega

    # Store values
    r_opt_list.append(r_opt)
    rpm_list.append(RPM_opt)
    pitch_list.append(pitch_opt)
    torque_list.append(torque_opt)
    power_list.append(P_opt)

    # Calculate RPM/Torque ratio and store it
    rpm_torque_ratio = RPM_opt / torque_opt
    rpm_torque_ratio_list.append(rpm_torque_ratio)

# Find the index where the RPM/Torque ratio is lowest
min_ratio_index = np.argmin(rpm_torque_ratio_list)
r_min_ratio = r_max_vals[min_ratio_index]
min_ratio_value = rpm_torque_ratio_list[min_ratio_index]
RPM_at_min_ratio = rpm_list[min_ratio_index]
torque_at_min_ratio = torque_list[min_ratio_index]

# Print the results for the lowest ratio
print(f"Minimum RPM/Torque ratio occurs at radius: {r_min_ratio} m with a ratio of {min_ratio_value}")
print(f"At this radius, the optimized RPM is: {RPM_at_min_ratio} RPM")
print(f"At this radius, the optimized Torque is: {torque_at_min_ratio} Nm")

# Plotting without Savitzky-Golay filter
plt.figure(figsize=(12, 8))

plt.subplot(3, 1, 1)
plt.plot(r_max_vals, rpm_list, label="RPM", color='blue')
plt.ylabel("Optimized RPM")
plt.grid(True)

plt.subplot(3, 1, 2)
plt.plot(r_max_vals, torque_list, label="Torque", color='green')
plt.ylabel("Torque [Nm]")
plt.grid(True)

plt.subplot(3, 1, 3)
plt.plot(r_max_vals, power_list, label="Power", color='red')
plt.xlabel("Upper Bound on Propeller Radius [m]")
plt.ylabel("Power [W]")
plt.grid(True)

plt.suptitle("Optimization of r and RPM vs Max Allowed Radius")
plt.tight_layout()
plt.show()
