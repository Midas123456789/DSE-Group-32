from ISA_Calculator import *
import numpy as np
import matplotlib.pyplot as plt

class WP_WS_Diagram:
    
    def __init__(self, V_cruise=0, V_stall=0, S=0, A=0, e=0, CD0=0, h=0, L=0, W=0, CL_MAX_clean=0, CL_MAX_landing=0, CL_MAX_TO=0, TOP=0, takeoff_distance=0):
        # Initialize parameters
        self.CL_max_clean = CL_MAX_clean
        self.CL_max_land  = CL_MAX_landing
        self.CL_max_TO    = CL_MAX_TO
        
        self.V_cruise = V_cruise
        self.V_stall = V_stall
        self.S = S
        self.A = A 
        self.e = e 
        self.CD0 = CD0
        self.h = h
        self.c = L
        self.W = W
        self.isa = ISA_Calculator(self.h, self.V_cruise, self.c)
        self.rho0 = 1.225
        self.rho = self.isa.results[self.h]["Density [kg/m³]"]
        self.TOP = TOP  # Takeoff Parameter (depends on aircraft class)
        self.takeoff_distance = takeoff_distance  # in meters
    

    
    def wing_loading(self, V, CL_max):
        """Calculate Wing Loading for given velocity and maximum lift coefficient"""
        return 0.5 * self.rho * V**2 * CL_max
    
    def take_off_loading(self, CL_TO_list=None):
        """
        Calculate W/P vs W/S curves for an array of CL_TO values.
        Returns a list of (W_P_range, W_S_range) tuples.
        """
        if CL_TO_list is None:
            CL_TO_list = [self.CL_max_TO / (1.1 ** 2)]

        W_S_range = np.linspace(1, 1500, 1000)
        results = []

        for CL_TO in CL_TO_list:
            effective_CL = CL_TO / (1.1 ** 2)
            W_P_range = (self.TOP / W_S_range) * effective_CL * (self.rho / self.rho0)
            results.append((W_P_range, W_S_range, CL_TO))
        
        return results
    
    def plot_wing_loading_constraints(self):
        """Plot the Wing Loading Constraints for Stall in Landing and Cruise"""
        
        x_max = 1000
        y_max = 1
        
        # Calculate wing loading for stall in landing condition
        WL_stall = self.wing_loading(self.V_stall, self.CL_max_land)
        
        # Calculate wing loading for cruise condition
        WL_cruise = self.wing_loading(self.V_cruise, self.CL_max_clean)
        
        # Calculate W/P vs W/S
        W_P_takeoff, W_S_takeoff = self.take_off_loading()
        
        # Plotting the constraints and W/P vs W/S
        plt.figure(figsize=(10, 6))
        
        # Plot Wing Loading Constraints
        plt.axvline(x=WL_stall, color='red', linestyle='--', label=f"Stall in Landing: {WL_stall:.2f} N/m²")
        plt.axvline(x=WL_cruise, color='blue', linestyle='--', label=f"Cruise: {WL_cruise:.2f} N/m²")
        
        # Fill the area left of the stall graph and above the takeoff curve in red
        plt.fill_between(W_S_takeoff, W_P_takeoff, y_max, color='red', alpha=0.3, label="Unfeasible Area")
        
        # Plot the W/P vs W/S relation
        plt.plot(W_S_takeoff, W_P_takeoff, label="W/P vs W/S (Takeoff)", color='green')
        
        plt.xlim(0, x_max)
        plt.ylim(0, y_max)
        plt.xlabel("Wing Loading (N/m²) / W/S (N/m²)")
        plt.ylabel("Power-to-Weight Ratio (W/P) [N/W]")
        plt.title("Wing Loading and Takeoff Power-to-Weight Relation")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

# Example Usage
if __name__ == "__main__":
    aircraft = WP_WS_Diagram(
        V_cruise       = 30,            # Cruise Speed (m/s)
        V_stall        = 10,            # Stall Speed (m/s)
        S              = 16,            # Wing Area (m²)
        A              = 7.5,           # Aspect Ratio
        e              = 0.80,          # Oswald Efficiency Factor
        CD0            = 0.030,         # Zero-lift Drag Coefficient
        h              = 0,         # Altitude (m)
        L              = 2,             # Wing Chord (m)
        W              = 6000,           # Weight (N)
        CL_MAX_clean   = 1.2,      # Clean Configuration Maximum Lift Coefficient
        CL_MAX_landing = 2.2,     # Landing Configuration Maximum Lift Coefficient
        CL_MAX_TO      = 2.0,       # Takeoff Configuration Maximum Lift Coefficient
        TOP            = 150,       # Takeoff Parameter
        takeoff_distance=500      # Takeoff Distance (meters)
    )
    
    aircraft.plot_wing_loading_constraints()
