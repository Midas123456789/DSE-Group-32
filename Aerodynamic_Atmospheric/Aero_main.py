# main.py
from Aircraft_Aerodynamics import AircraftAerodynamic
from Aero_plotting import plot_feasible_S_V, plot_A_LD, plot_h_V, plot_P_V, plot_h_V, plot_h_W
from Aero_optim import optimize_velocity_for_min_drag, optimize_altitude_for_min_drag, parameter_sweep, determine_states_altitude, calculate_W, calculate_W_bat
from Power_calculation import sweep_propeller_radius
import numpy as np
from math import inf
import aerosandbox as asb
import matplotlib.pyplot as plt

if __name__ == "__main__":




    


    plot_P_V()
    plot_h_V()
    plot_h_W()

    # states = determine_states_altitude(aircraft, np.arange(0, 30000, 100))
    # for i, state in enumerate(states):
    #     h, rho, V, D, L, M = state
    #     opt_state = sweep_propeller_radius(asb=asb, D=D, rho=rho, v_o=V, r_min=0.75, r_max=6, r_step=0.05)
    #     states[i] = [h, rho, V, D, L, M, opt_state[0], opt_state[1]]
    # states = np.array(states)



    # # Plot Density vs Power
    # plt.figure()
    # plt.plot(states[:, 0], states[:, 7], label='Density vs Power')
    # plt.xlabel('Altitude (m)')
    # plt.ylabel('Power (W)')
    # plt.title('Altitude vs Power')
    # plt.legend()
    # plt.grid()

    # # Plot Altitude vs Velocity
    # plt.figure()
    # plt.plot(states[:, 0], states[:, 2], label='Altitude vs Velocity')
    # plt.xlabel('Altitude (m)')
    # plt.ylabel('Velocity (m/s)')
    # plt.title('Altitude vs Velocity')
    # plt.legend()
    # plt.grid()

    # # Plot Altitude vs Drag
    # plt.figure()
    # plt.plot(states[:, 0], states[:, 3], label='Altitude vs Drag')
    # plt.xlabel('Altitude (m)')
    # plt.ylabel('Drag (N)')
    # plt.title('Altitude vs Drag')
    # plt.legend()
    # plt.grid()

    # # Plot Altitude vs Lift
    # plt.figure()
    # plt.plot(states[:, 0], states[:, 4], label='Altitude vs Lift')
    # plt.xlabel('Altitude (m)')
    # plt.ylabel('Lift (N)')
    # plt.title('Altitude vs Lift')
    # plt.legend()
    # plt.grid()

    # # Plot Altitude vs Mach Number
    # plt.figure()
    # plt.plot(states[:, 0], states[:, 5], label='Altitude vs Mach Number')
    # plt.xlabel('Altitude (m)')
    # plt.ylabel('Mach Number')
    # plt.title('Altitude vs Mach Number')
    # plt.legend()
    # plt.grid()
    # plt.show()




    # # print(f"Optimal Altitude: {optimal_h:.2f} m, Drag: {drag:.2f} N")
