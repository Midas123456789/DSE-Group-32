import matplotlib.pyplot as plt
import numpy as np

def plot_feasible_S_V(aircraft, CL_list):
    V = np.linspace(10, 200, 500)
    plt.figure(figsize=(10, 6))
    for cl in CL_list:
        if cl < 0:
            raise ValueError(f"Lift coefficient must be > 0. Got {cl}.")
        S = (2 * aircraft.weight) / (aircraft.rho * V**2 * cl)
        plt.plot(V, S, label=f'CL = {cl}')
        plt.fill_between(V, 0, S, alpha=0.1)
    plt.axhline(30, color='red', linestyle='--', label='30 m² Wing Area')
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Wing Area (m²)')
    plt.title(f'Feasible Wing Area vs Velocity\nAltitude: {aircraft.altitude} m | Weight: {aircraft.weight} N')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_h_V(aircraft):
    altitude_range = aircraft.altitude_range
    V = []
    for h in altitude_range:
        rho = aircraft.altitude_data[h]["Density [kg/m³]"]
        V.append(np.sqrt((2 * aircraft.weight) / (rho * aircraft.S * aircraft.CL)))
    plt.figure(figsize=(10, 6))
    plt.plot(V, altitude_range, label='Max Altitude')
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Altitude (m)')
    plt.title(f'Altitude vs Velocity\nWeight: {aircraft.weight} N')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_A_LD(aircraft, CL_list):
    A_range = np.linspace(1, 20, 500)
    plt.figure(figsize=(10, 6))
    for cl in CL_list:
        CD = aircraft.drag_polar(CL=cl, A=A_range)
        LD = aircraft.Lift(CL=cl) / aircraft.Drag(CD=CD)
        plt.plot(A_range, LD, label=f'CL = {cl}')
    plt.xlabel('Aspect Ratio (A)')
    plt.ylabel('Lift-to-Drag Ratio (L/D)')
    plt.title('L/D Ratio vs Aspect Ratio')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()