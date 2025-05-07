import matplotlib.pyplot as plt
import numpy as np
from Aerodynamics import *

class AircraftAerodynamic(Aerodynamic):
    
    def __init__(self, W=0, h=0, V=0, S=0, A=0, e=0, CD0=0, CL=0):
        super().__init__(weight=W, altitude=h)
        
        self.rho = self.altitude_data[self.altitude]["Density [kg/m³]"]
        self.V = V
        self.S = S
        self.A = A 
        self.e = e 
        self.CD0 = CD0
        self.CL = CL
    
    def drag_polar(self):
        self.CD = self.CD0 + ((self.CL)**2)/(np.pi * self.A * self.e)
        self.D = (1/2) * self.rho * (self.V ** 2) * self.S
        return self.D
    
    def plot_feasible_S_V(self, CL_list):
        """
        Plot feasible wing area (S) vs. velocity (V) for different coefficients of lift.
        
        Parameters:
        - CL: float or list of floats (lift coefficients)
        """
        
        # Define velocity range
        V = np.linspace(10, 200, 500)
        
        # Normalize CL input
        if not isinstance(CL_list, (list, tuple, np.ndarray)):
            CL_list = [CL_list]
        
        # Create plot
        plt.figure(figsize=(10, 6))
        
        for cl in CL_list:
            if cl < 0:
                raise ValueError(f"Lift coefficient must be bigger than 0. Got {cl}.")
            S = (2 * self.weight) / (self.rho * V**2 * cl)
            plt.plot(V, S, label=f'CL = {cl}')
            plt.fill_between(V, 0, S, alpha=0.1)
        
        # Add labels and formatting
        plt.xlabel('Velocity (m/s)')
        plt.ylabel('Wing Area (m²)')
        plt.title(f'Feasible Wing Area vs Velocity\nAltitude: {self.altitude} m | Weight: {self.weight} N')
        plt.xlim(10, 200)
        plt.ylim(bottom=0)
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()


# Example Usage
if __name__ == "__main__":
    aircraft = AircraftAerodynamic(
        W=4000, 
        h=20000, 
        V=25, 
        S=10, 
        A=12, 
        e=0.9, 
        CD0=0,
        CL=2
    )
    
    aircraft.plot_feasible_S_V(CL_list=[0.5, 1, 1.5, 2, 2.5, 3])
    print(aircraft.drag_polar())