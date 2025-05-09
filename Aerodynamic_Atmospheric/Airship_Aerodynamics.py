import numpy as np
import matplotlib.pyplot as plt
from Aerodynamics import Aerodynamic

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
        
        molarmass_gas = {
            "hydrogen": 0.002016,  # kg/mol
            "helium": 0.0040026    # kg/mol
        }[self.gas_type]


        V = []
        for i in self.altitude_range:
            rho_gas = (self.altitude_data[i]["Pressure [Pa]"] * molarmass_gas) / (self.R * self.altitude_data[i]["Temperature [K]"])
            V.append(self.weight / (self.altitude_data[i]["Density [kg/m³]"] - rho_gas))  # Convert to m³

        plt.figure(figsize=(10, 6))
        plt.plot(self.altitude_range, V, label=f'Gas Type: {self.gas_type.capitalize()}')
        plt.xlabel('Altitude (m)')
        plt.ylabel('Volume (m³)')
        plt.title(f'Feasible Altitude (h) vs Volume (V)\nWeight: {self.weight} N')
        plt.xlim(0, 25000)
        #plt.axhline(30000, color='red', linestyle='--', label='Max Volume (1e6 ft³)')
        plt.ylim(0, 100000)
        plt.grid(True)
        plt.legend()
        plt.show()


# Airship example
airship = AirshipAerodynamic(weight=9000, altitude=20000, gas_type="hydrogen")
airship.plot_feasible_h_V()
