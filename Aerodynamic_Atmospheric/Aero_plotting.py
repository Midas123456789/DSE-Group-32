import matplotlib.pyplot as plt
import numpy as np
import aerosandbox as asb
from Power_calculation import sweep_propeller_radius
from Aircraft_Aerodynamics import AircraftAerodynamic

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

def plot_P_V():
    Vrange = np.arange(10, 100, 1)

    P_lst = []
    for V in Vrange:
        aircraft = AircraftAerodynamic(W=4000, h=15000, V=V, S=30, A=25, e=0.7, CD0=0.04, CL=1.2
        )
        rho = aircraft.rho
        V = aircraft.V
        D = aircraft.Drag()

        P = sweep_propeller_radius(asb=asb, D=D, rho=rho, v_o=V)[1]

        if V == 20 or V == 50:
           specialpower = P
           otherpower = np.sqrt(2*aircraft.weight**3*aircraft.CD**2/(aircraft.S*aircraft.rho*aircraft.CL**3))
           Drag = D
           Mach = aircraft.M
           CD = aircraft.CD
        P_lst.append(P)



    plt.figure(figsize=(10, 6))
    plt.plot(Vrange, np.array(P_lst)/1000, label='Power vs Velocity')
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Power (kW)')
    # plt.title(f'Power vs Velocity\nAltitude: {h} m | Weight: {aircraft.weight} N')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
    print('hello')
    print(specialpower)
    print(otherpower)
    print(Drag)
    print(Mach)
    print(CD)
    print(Vrange)




