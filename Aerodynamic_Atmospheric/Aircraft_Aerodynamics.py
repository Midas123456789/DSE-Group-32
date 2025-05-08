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
        self.CD = self.drag_polar()
    
    def drag_polar(self, CL=None, CD0=None, A=None, e=None):
        """
        Calculate the drag polar (CD) based on the given or default parameters.

        Parameters:
        - CL: Coefficient of lift (optional, overrides self.CL if provided)
        - CD0: Zero-lift drag coefficient (optional, overrides self.CD0 if provided)
        - A: Aspect ratio (optional, overrides self.A if provided)
        - e: Oswald efficiency factor (optional, overrides self.e if provided)

        Returns:
        - CD: Coefficient of drag
        """
        CL = CL if CL is not None else self.CL
        CD0 = CD0 if CD0 is not None else self.CD0
        A = A if A is not None else self.A
        e = e if e is not None else self.e

        CD = CD0 + (CL ** 2) / (np.pi * A * e)
        return CD
    
    def Drag(self, CD=None, V=None, S=None, rho=None):
        """
        Calculate drag force based on the drag equation.
        
        Parameters:
        - CD: Coefficient of drag
        - V: Velocity in m/s
        - S: Reference area in m²
        
        Returns:
        - Drag force in Newtons
        """

        CD = CD if CD is not None else self.CD
        V = V if V is not None else self.V
        S = S if S is not None else self.S
        rho = rho if rho is not None else self.rho

        
        return 0.5 * rho * (V ** 2) * S * CD
    
    def Lift(self, CL=0, V=0, S=0, rho=0):
        """
        Calculate lift force based on the lift equation.
        
        Parameters:
        - CL: Coefficient of lift
        - V: Velocity in m/s
        - S: Wing area in m²
        
        Returns:
        - Lift force in Newtons
        """
        CL = CL if CL is not None else self.CL
        V = V if V is not None else self.V
        S = S if S is not None else self.S
        rho = rho if rho is not None else self.rho

        return 0.5 * rho * (V ** 2) * S * CL
    
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
        plt.axhline(30, color='red', linestyle='--', label='30 m² Wing Area')
        plt.title(f'Feasible Wing Area vs Velocity\nAltitude: {self.altitude} m | Weight: {self.weight} N')
        plt.xlim(10, 200)
        plt.ylim(0, 100)
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()
    
    def plot_h_V(self):
        """
        Plot feasible altitude (h) vs. velocity (V) for the aircraft.
        
        Parameters:
        - None
        """
        
        # Define altitude range
        altitude_range = self.altitude_range
        V = []

        for h in altitude_range:
            rho = self.altitude_data[h]["Density [kg/m³]"]
            V.append(np.sqrt((2 * self.weight) / (rho * self.S * self.CL)))
        
        # Create plot
        plt.figure(figsize=(10, 6))
        plt.plot(V, altitude_range, label='Maximum Altitude')
        
        # Add labels and formatting
        plt.xlabel('Velocity (m/s)')
        plt.ylabel('Altitude (m)')
        plt.title(f'Feasible Altitude vs Velocity\nWeight: {self.weight} N')
        plt.xlim(10, 200)
        plt.ylim(bottom=0)
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()

    def plot_A_LD(self, CL_list):
        """
        Plot feasible wing area (S) vs. velocity (V) for different coefficients of lift.

        Parameters:
        - CL: float or list of floats (lift coefficients)
        """


        if not isinstance(CL_list, list):
            CL_list = [CL_list]
        
        # Define Aspect Ratio (A) range
        A_range = np.linspace(1, 20, 500)

        CD = self.drag_polar(CL=self.CL, A=A_range, e=self.e, CD0=self.CD0)
        LD = self.Lift(CL=self.CL, V=self.V, S=self.S, rho=self.rho) / self.Drag(CD=CD, V=self.V, S=self.S, rho=self.rho)

        plt.plot(A_range, LD, label=f'CL = {self.CL}')
        plt.xlabel('Aspect Ratio (A)')
        plt.ylabel('Lift-to-Drag Ratio (L/D)')
        plt.title(f'Lift-to-Drag Ratio vs Aspect Ratio\nAltitude: {self.altitude} m | Weight: {self.weight} N')
        plt.xlim(1, 20)
        plt.ylim(bottom=0)
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()


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


# Example Usage
if __name__ == "__main__":
    aircraft = AircraftAerodynamic(
        W=5000, 
        h=20000, 
        V=50, 
        S=30, 
        A=12, 
        e=0.75, 
        CD0=0.05,
        CL=2
    )
    
    aircraft.plot_feasible_S_V(CL_list=[2, 2.5, 3])
    aircraft.plot_A_LD()
    aircraft.plot_h_V()
    CD = aircraft.drag_polar()
    drag = aircraft.Drag(CD)
    print(f"Drag: {drag:.2f} N")