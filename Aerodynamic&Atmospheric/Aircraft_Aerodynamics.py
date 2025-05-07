import numpy as np
import matplotlib.pyplot as plt
from Aerodynamics import Aerodynamic

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





    def plot_feasible_A_LD(self, CL):
        """
        Plot feasible wing area (S) vs. velocity (V) for different coefficients of lift.

        Parameters:
        - CL: float or list of floats (lift coefficients)
        """


        if not isinstance(CL, list):
            CL = [CL]

        # input A, e, CL, CD0, S

        A = np.linspace(10, 100, 500)

    def WP_WS_Diagram(self, Vs= 10, CL_max = 2, CL_TO = 2.4, TOP = 500, sigma = 1):
        """
        Wing loading diagram
        """

        WS_range = np.linspace(1, 10000, 500)

        WSmax = 0.5 * self.altitude_data[self.altitude]["Density [kg/m³]"] * Vs**2 * CL_max

        PW_TO = TOP/(WS_range * CL_TO * sigma)

        plt.figure(figsize=(10, 6))
        plt.axvline(WSmax, color='red', linestyle='--', label='Max Wing Loading')
        plt.plot(WS_range, PW_TO, label='Power-to-Weight Ratio (P/W)')
        plt.fill_between(WS_range, 0, PW_TO, alpha=0.1)
        plt.xlabel('Wing Loading (W/S)')
        plt.xlim(0, 100)
        plt.ylabel('(P/W)')
        plt.show()





# Aircraft example
aircraft = AircraftAerodynamic(weight=4000, altitude=20000)
# aircraft.plot_feasible_S_V([1, 1.5, 2, 2.5, 3])
aircraft.WP_WS_Diagram()