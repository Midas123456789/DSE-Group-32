import matplotlib.pyplot as plt
import numpy as np
from Power_calculation import compute_propeller_power

def plot_feasible_S_V(aircraft, CL_list):
    V = np.linspace(10, 200, 500)
    plt.figure(figsize=(10, 6))
    for cl in CL_list:
        if cl <= 0:
            raise ValueError(f"Lift coefficient must be > 0. Got {cl}.")
        # Use inputs from the new aircraft class
        S = (2 * aircraft.class_I.estimated_MTOW) / (aircraft.aero.altitude_data[aircraft.inputs.h_cruise]["Density [kg/m3]"] * V**2 * cl)
        plt.plot(V, S, label=f'CL = {cl}')
        plt.fill_between(V, 0, S, alpha=0.1)
    plt.axhline(aircraft.inputs.S, color='red', linestyle='--', label=f'{aircraft.inputs.S} m² Wing Area')
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Wing Area (m²)')
    plt.title(f'Feasible Wing Area vs Velocity\nAltitude: {aircraft.inputs.h_cruise} m | Weight: {aircraft.class_I.estimated_MTOW:.0f} N')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_h_V(aircraft):
    altitude_range = np.arange(0, 30000, 100)
    Vmin = []
    Vmax = []
    Pmax = 2000  # Watts

    for h in altitude_range:
        rho = aircraft.aero.altitude_data[h]["Density [kg/m3]"]
        Vmin.append(np.sqrt((2 * aircraft.class_I.estimated_MTOW) / (rho * aircraft.inputs.S * aircraft.inputs.CL_cruise)))

        CD = aircraft.aero.drag_polar()
        Vmax.append((Pmax * 2 / (rho * aircraft.inputs.S * CD)) ** (1 / 3))

    plt.figure(figsize=(10, 6))
    plt.plot(Vmin, altitude_range, label='Min Speed')
    plt.plot(Vmax, altitude_range, label='Max Speed')
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Altitude (m)')
    plt.title(f'Altitude vs Velocity\nWeight: {aircraft.class_I.estimated_MTOW:.0f} N')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_A_LD(aircraft, CL_list):
    A_range = np.linspace(1, 20, 500)
    plt.figure(figsize=(10, 6))
    for cl in CL_list:
        # Use drag polar method from the new aircraft class
        CD = aircraft.aero.drag_polar(CL=cl, A=A_range)
        LD = aircraft.aero.Lift(CL=cl) / aircraft.aero.Drag(CD=CD)
        plt.plot(A_range, LD, label=f'CL = {cl}')
    plt.xlabel('Aspect Ratio (A)')
    plt.ylabel('Lift-to-Drag Ratio (L/D)')
    plt.title('L/D Ratio vs Aspect Ratio')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_P_V(aircraft):
    Vrange = np.arange(10, 60, 1)
    hrange = np.arange(15000, 25000, 5000)

    plt.figure(figsize=(10, 6))
    for h in hrange:
        P_lst = []
        P2_lst = []
        rho = aircraft.aero.altitude_data[h]["Density [kg/m3]"]
        for V in Vrange:
            # Use drag method and other aerodynamic properties from the new class
            D = aircraft.aero.Drag(V=V, rho=rho)
            print(D)
            r = 2  # Example prop radius
            P = compute_propeller_power(D, rho, V, r)
            P2 = 0.5 * rho * V**3 * aircraft.inputs.S * aircraft.aero.CD
            P_lst.append(P)
            P2_lst.append(P2)

        plt.plot(Vrange, np.array(P_lst) / 1000, label=f'Rotor Power at {h} m')
        plt.plot(Vrange, np.array(P2_lst) / 1000, '--', label=f'Aerodynamic Power at {h} m')

    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Power (kW)')
    plt.title('Power vs Velocity at Different Altitudes')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_h_W(aircraft):
    # This function seems to be unrelated to the `AircraftAerodynamic` class
    # I would assume you'd want to keep the legacy code for this one
    h_list = np.arange(0, 20000, 100)
    W_lst = []

    for h in h_list:
        rho = aircraft.aero.altitude_data[h]["Density [kg/m3]"]
        weight = 0
        for _ in range(100):
            V = aircraft.aero.min_velocity(rho=rho, weight=weight)
            D = aircraft.aero.Drag(rho=rho, V=V)
            Preq = D * V
            Wbat = calculate_W_bat(Preq)
            weight = calculate_W(Wbat)
        if weight != float('inf'):
            W_lst.append(weight)
        else:
            W_lst.append(0)

    plt.plot(h_list, W_lst, label='Weight vs Altitude')
    plt.xlabel('Altitude (m)')
    plt.ylabel('Weight (N)')
    plt.ylim(0, 20000)
    plt.title('Weight vs Altitude')
    plt.legend()
    plt.grid()
    plt.show()

def calculate_W_bat(P):
    # Assume 18kW payload, 8 hour night, 450Wh/kg battery
    return ((P + 18000) * 8 / 435) * 9.81

def calculate_W(Wbat):
    # Assume 100 kg payload, 1/10 payload structure weight
    return Wbat + 10 * 100


def plot_h_Preq(aircraft):
    h_list = np.arange(0, 12000, 100)
    Preq_lst = []
    V_lst = []

    for h in h_list:
        rho = aircraft.aero.altitude_data[h]["Density [kg/m3]"]
        weight = aircraft.class_I.estimated_MTOW
        for _ in range(100):
            V = aircraft.aero.min_velocity(rho=rho, weight=weight)
            D = aircraft.aero.Drag(rho=rho, V=V)
            Preq = D * V
            Wbat = calculate_W_bat(Preq)
            weight = calculate_W(Wbat)
        if weight != float('inf'):
            Preq_lst.append(Preq)
            V_lst.append(V)


    fig, ax1 = plt.subplots()

    ax1.plot(h_list, Preq_lst, label='Power vs Altitude', color='blue')
    ax1.set_xlabel('Altitude (m)')
    ax1.set_ylabel('Power (W)', color='blue')
    ax1.tick_params(axis='y', labelcolor='blue')

    ax2 = ax1.twinx()
    ax2.plot(h_list, V_lst, label='Velocity vs Altitude', linestyle='--', color='green')
    ax2.set_ylabel('Velocity (m/s)', color='green')
    ax2.tick_params(axis='y', labelcolor='green')

    fig.tight_layout()
    plt.title('Power and Velocity vs Altitude')
    plt.xlabel('Altitude (m)')
    plt.ylabel('Power (W)')
    plt.title('Power vs Altitude')
    plt.legend()
    plt.grid()
    plt.show()


# Local testing
if __name__ == "__main__":
    pass