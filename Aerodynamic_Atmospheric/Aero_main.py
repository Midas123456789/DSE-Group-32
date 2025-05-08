# main.py
from Aircraft_Aerodynamics import AircraftAerodynamic
from Aero_plotting import plot_feasible_S_V, plot_A_LD, plot_h_V
from Aero_optim import optimize_velocity_for_min_drag, optimize_altitude_for_min_drag, parameter_sweep, determine_states_altitude
from Power_calculation import sweep_propeller_radius
import numpy as np
import aerosandbox as asb
import matplotlib.pyplot as plt

if __name__ == "__main__":
    aircraft = AircraftAerodynamic(W=7000, h=20000, V=20, S=50, A=12, e=1, CD0=0.1, CL=3)

    states = determine_states_altitude(aircraft, np.arange(0, 30000, 100))
    for i, state in enumerate(states):
        h, rho, V, D, L, M = state
        # opt_state = sweep_propeller_radius(asb=asb, D=D, rho=rho, v_o=V, r_min=0.75, r_max=6, r_step=0.05)
        # states[i] = [h, rho, V, D, L, opt_state[0], opt_state[1]]
    states = np.array(states)
    print(states)

    plt.plot(states[:, 0], states[:, 3], label='Density vs Altitude')
    plt.show()




    # print(f"Optimal Altitude: {optimal_h:.2f} m, Drag: {drag:.2f} N")
