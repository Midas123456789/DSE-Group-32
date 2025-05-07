import numpy as np
import matplotlib.pyplot as plt
from ISA_Calculator import ISA_Calculator

class Aerodynamic:
    def __init__(self, weight, altitude):
        """
        Base aerodynamic model.
        Parameters:
        - weight: in Newtons
        - altitude: in meters
        """
        self.weight = weight
        self.altitude = altitude

        self.altitude_range = np.linspace(0, 30000, 100).tolist()
        self.isa = ISA_Calculator(self.altitude_range)
        self.altitude_data = self.isa.results


class AircraftAerodynamic(Aerodynamic):
    def plot_feasible_S_V(self, CL):
        """
        Plot feasible wing area (S) vs. velocity (V) for different coefficients of lift.

        Parameters:
        - CL: float or list of floats (lift coefficients)
        """
        rho = self.altitude_data[self.altitude]["Density [kg/m³]"]
        V = np.linspace(10, 200, 500)

        if not isinstance(CL, list):
            CL = [CL]

        plt.figure(figsize=(10, 6))
        for cl in CL:
            S = (2 * self.weight) / (rho * V**2 * cl)
            plt.plot(V, S, label=f'CL = {cl}')
            plt.fill_between(V, 0, S, alpha=0.1)

        plt.xlabel('Velocity (m/s)')
        plt.ylabel('Wing Area (m²)')
        plt.title(f'Feasible Wing Area (S) vs Velocity (V)\nAltitude: {self.altitude} m, Weight: {self.weight} N')
        plt.xlim(0, 100)
        plt.grid(True)
        plt.legend()
        plt.show()

    def estimate_drag(self, S, V, CD0, k):
        """
        Estimate drag using the drag equation.

        Parameters:
        - S: wing area in m²
        - V: velocity in m/s
        - CD0: zero-lift drag coefficient
        - k: induced drag factor

        Returns:
        - Drag force in Newtons
        """
        rho = self.rho
        q = 0.5 * rho * V**2
        CD = CD0 + k * (self.weight / (0.5 * rho * V**2 * S))**2

class AirshipAerodynamic(Aerodynamic):
    def __init__(self, weight, altitude, gas_type):
        """
        Airship aerodynamic model.

        Parameters:
        - gas_type: 'hydrogen' or 'helium'
        """
        super().__init__(weight, altitude)

        if gas_type not in ["hydrogen", "helium"]:
            raise ValueError("Gas type must be 'hydrogen' or 'helium'")
        self.gas_type = gas_type

    def plot_feasible_h_V(self):
        """
        Plot feasible altitude (h) vs. Volume for the airship.
        L>W

        Parameters:
        - None
        """
        
        rho_gas = {
            "hydrogen": 0.0899,  # kg/m³
            "helium": 0.1786     # kg/m³
        }[self.gas_type]

        V = []
        for i in self.altitude_range:
            V.append(self.weight / (self.altitude_data[i]["Density [kg/m³]"] - rho_gas))  # Convert to m³

        plt.figure(figsize=(10, 6))
        plt.plot(self.altitude_range, V, label=f'Gas Type: {self.gas_type.capitalize()}')
        # plt.fill_between(h, 0, V, alpha=0.1)
        plt.xlabel('Altitude (m)')
        plt.ylabel('Volume (m³)')
        plt.title(f'Feasible Altitude (h) vs Volume (V)\nWeight: {self.weight} N')
        plt.xlim(0, 25000)
        plt.axhline(30000, color='red', linestyle='--', label='Max Volume (1e6 ft³)')
        plt.ylim(0, 35000)
        plt.grid(True)
        plt.legend()
        plt.show()

# Aircraft example
aircraft = AircraftAerodynamic(weight=4000, altitude=20000)
aircraft.plot_feasible_S_V([1, 1.5, 2, 2.5, 3])

# Airship example
airship = AirshipAerodynamic(weight=4000, altitude=20000, gas_type="hydrogen")
airship.plot_feasible_h_V()
# print("Buoyant Force:", airship.bouyant_force(V=500))
