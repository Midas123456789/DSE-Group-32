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
    hrange = np.arange(0, 30000, 5000)

    for h in hrange:
        P_lst = []
        P2_lst = []
        P3_lst = []
        for V in Vrange:
            aircraft = AircraftAerodynamic(W=4000, h=h, V=V, S=30, A=25, e=0.7, CD0=0.04, CL=1.2
            )
            rho = aircraft.rho
            V = aircraft.V
            D = aircraft.Drag()

            P = sweep_propeller_radius(asb=asb, D=D, rho=rho, v_o=V)[1]
            P2 = 0.5 * rho * (V ** 3) * aircraft.S * aircraft.CD

            P_lst.append(P)
            P2_lst.append(P2)
        
        plt.plot(Vrange, np.array(P_lst)/1000, label=f'Power vs Velocity (Rotor) at {h} m')
        plt.plot(Vrange, np.array(P2_lst)/1000, label=f'Power vs Velocity (Aerodynamic) at {h} m', linestyle='--')





    plt.figure(figsize=(10, 6))
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Power (kW)')
    # plt.title(f'Power vs Velocity\nAltitude: {h} m | Weight: {aircraft.weight} N')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_h_V():
    altitude_range = np.arange(0, 30000, 100)
    Vmin = []
    for h in altitude_range:
        aircraft = AircraftAerodynamic(W=6000, h=h, V=0, S=30, A=25, e=0.7, CD0=0.04, CL=1.2)
        rho = aircraft.altitude_data[h]["Density [kg/m³]"]
        Vmin.append(np.sqrt((2 * aircraft.weight) / (rho * aircraft.S * aircraft.CL)))

    Vmax = []
    Pmax = 2000 # Watts
    for h in altitude_range:
        aircraft = AircraftAerodynamic(W=6000, h=h, V=0, S=30, A=25, e=0.7, CD0=0.04, CL=1.2)
        rho = aircraft.altitude_data[h]["Density [kg/m³]"]
        CD = aircraft.drag_polar()
        Vmax.append((Pmax*2/(rho*aircraft.S*CD))**(1/3))
    
    plt.figure(figsize=(10, 6))
    plt.plot(Vmin, altitude_range, label='Min Speed')
    plt.plot(Vmax, altitude_range, label='Max Speed')
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Altitude (m)')
    plt.title(f'Altitude vs Velocity\nWeight: {aircraft.weight} N')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_h_W():
    h_list = np.arange(0, 30000, 100)
    W_lst = []
    V_lst = []

    for h in h_list:
        aircraft = AircraftAerodynamic(W=4000, h=h, V=0, S=30, A=20, e=0.7, CD0=0.04, CL=1.2)

        for i in range(0, 100):
            aircraft.V = aircraft.min_velocity()
            D = aircraft.Drag()
            Preq = D * aircraft.V
            Wbat = calculate_W_bat(Preq)
            aircraft.weight = calculate_W(Wbat)
        if aircraft.weight != float('inf'):
            W_lst.append(aircraft.weight)
            V_lst.append(aircraft.V)
        else:
            W_lst.append(0)
            V_lst.append(0)

    plt.plot(h_list, W_lst, label='Weight vs Altitude')
    plt.xlabel('Altitude (m)')
    plt.ylabel('Weight (N)')
    plt.ylim(0, 20000)
    plt.title('Weight vs Altitude')
    plt.legend()
    plt.grid()
    plt.show()

def calculate_W_bat(P):
# Assume 18kW payload
# Assume 8 hour night
# Assume 450Wh/kg battery
    return ((P + 18000) * 8 / 435)*9.81

def calculate_W(Wbat):
    # Assume 100 kg payload
    # Assume 1/10 payload structure weight
    return Wbat + 10*100
