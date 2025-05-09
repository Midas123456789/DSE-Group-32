import aerosandbox as asb
import aerosandbox.numpy as np
from aerosandbox import cas
from numpy import pi
from aerosandbox import Atmosphere
import matplotlib.pyplot as plt
from Aircraft_Aerodynamics import AircraftAerodynamic

# Constants
# altitude = 20000  # meters
# D = 400           # Thrust [N]
# v_o = 40          # Flight speed [m/s]
# atmo = Atmosphere(altitude)
# rho = atmo.density()

def sweep_propeller_radius(asb, D, rho, v_o, r_min=0.75, r_max=6, r_step=0.05):
    """
    Sweep the upper bound of the propeller radius and optimize for minimum power at each step.
    
    Parameters:
    - asb: module or object with `Opti()` method (e.g., from aerosandbox)
    - D: thrust (N)
    - rho: air density (kg/m^3)
    - v_o: freestream velocity (m/s)
    - r_min, r_max, r_step: sweep range for upper bound on propeller radius (in meters)

    Returns:
    - Dictionary with lists of optimized values and plot is shown.
    """
    
    r_max_vals = np.arange(r_min, r_max, r_step)

    # Storage lists
    r_opt_list = []
    rpm_list = []
    pitch_list = []
    torque_list = []
    power_list = []
    rpm_torque_ratio_list = []

    opti = asb.Opti()
    r = 2

    Adisk = pi * r**2
    T = D
    power = 0.5 * T * v_o * (cas.sqrt(T / (0.5 * rho * Adisk * v_o**2) + 1) + 1)

    opti.minimize(power)
    sol = opti.solve()
    r_opt = sol.value(r)
    P_opt = sol.value(power)



    
    # for r_max in r_max_vals:
    #     opti = asb.Opti()

    #     # Variables
    #     r = opti.variable(init_guess=1, lower_bound=0.5, upper_bound=r_max)
    #     RPM = opti.variable(init_guess=600, lower_bound=500, upper_bound=20000)

    #     # Derived quantities
    #     Adisk = pi * r**2
    #     T = D
    #     power = 0.5 * T * v_o * (cas.sqrt(T / (0.5 * rho * Adisk * v_o**2) + 1) + 1)
    #     pitch = v_o / (RPM / 60)

    #     # Objective
    #     opti.minimize(power)

    #     # Solve
    #     sol = opti.solve()
    #     r_opt = sol.value(r)
    #     RPM_opt = sol.value(RPM)
    #     P_opt = sol.value(power)
    #     pitch_opt = sol.value(pitch)
    #     omega = 2 * pi * RPM_opt / 60
    #     torque_opt = P_opt / omega

    #     # Store results
    #     r_opt_list.append(r_opt)
    #     rpm_list.append(RPM_opt)
    #     pitch_list.append(pitch_opt)
    #     torque_list.append(torque_opt)
    #     power_list.append(P_opt)
    #     rpm_torque_ratio_list.append(RPM_opt / torque_opt)

    # # Find lowest RPM/Torque ratio
    # min_ratio_index = np.argmin(rpm_torque_ratio_list)
    # r_min_ratio = r_max_vals[min_ratio_index]
    # min_ratio_value = rpm_torque_ratio_list[min_ratio_index]
    # RPM_at_min_ratio = rpm_list[min_ratio_index]
    # torque_at_min_ratio = torque_list[min_ratio_index]

    # print(f"Minimum RPM/Torque ratio occurs at radius: {r_min_ratio:.2f} m with a ratio of {min_ratio_value:.2f}")
    # print(f"At this radius, the optimized RPM is: {RPM_at_min_ratio:.2f} RPM")
    # print(f"At this radius, the optimized Torque is: {torque_at_min_ratio:.2f} Nm")

    # # Plot
    # plt.figure(figsize=(12, 8))

    # plt.subplot(3, 1, 1)
    # plt.plot(r_max_vals, rpm_list, label="RPM", color='blue')
    # plt.ylabel("Optimized RPM")
    # plt.grid(True)

    # plt.subplot(3, 1, 2)
    # plt.plot(r_max_vals, torque_list, label="Torque", color='green')
    # plt.ylabel("Torque [Nm]")
    # plt.grid(True)

    # plt.subplot(3, 1, 3)
    # plt.plot(r_max_vals, power_list, label="Power", color='red')
    # plt.xlabel("Upper Bound on Propeller Radius [m]")
    # plt.ylabel("Power [W]")
    # plt.grid(True)

    # plt.suptitle("Optimization of r and RPM vs Max Allowed Radius")
    # plt.tight_layout()
    # plt.show()
    return [r_opt, P_opt]

    return {
        "r_max_vals": r_max_vals,
        "r_opt_list": r_opt_list,
        "rpm_list": rpm_list,
        "pitch_list": pitch_list,
        "torque_list": torque_list,
        "power_list": power_list,
        "rpm_torque_ratio_list": rpm_torque_ratio_list,
        "min_ratio_info": {
            "radius": r_min_ratio,
            "ratio": min_ratio_value,
            "RPM": RPM_at_min_ratio,
            "torque": torque_at_min_ratio
        }
    }


aircraft = AircraftAerodynamic(W=6000, h=15000, V=20, S=30, A=25, e=0.7, CD0=0.04, CL=1.2
        )

D = aircraft.Drag()
rho = aircraft.rho
v_o = aircraft.V

results = sweep_propeller_radius(asb=asb, D=D, rho=rho, v_o=v_o)
print(D, rho, v_o)
print(results)