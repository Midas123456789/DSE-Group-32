from ISA_Calculator import *
import numpy as np
import matplotlib.pyplot as plt

class WP_WS_Diagram:
    
    def __init__(self, V_cruise=0, V_stall=0, S=0, A=0, e=0, CD0=0, h=0, L=0, W=0,
            CL_MAX_clean=None, CL_MAX_land=0, CL_MAX_TO=None, prop_setting=0, TOP=0, landing_distance=0, f=1, n_p=1):
        # Ensure clean and TO CLs are lists
        self.x_max = 1700
        self.y_max = 1
        
        self.CL_max_clean = CL_MAX_clean if isinstance(CL_MAX_clean, list) else [CL_MAX_clean]
        self.CL_max_land  = CL_MAX_land if isinstance(CL_MAX_land, list) else [CL_MAX_land]
        self.CL_max_TO    = CL_MAX_TO if isinstance(CL_MAX_TO, list) else [CL_MAX_TO]
        self.prop_setting = prop_setting if isinstance(prop_setting, list) else [prop_setting]
        
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
        self.landing_distance = landing_distance
        self.f = f
        self.n_p = n_p
    
    def wing_loading(self, V, CL_max):
        """Calculate Wing Loading for given velocity and maximum lift coefficient in cruise"""
        return 0.5 * self.rho * V**2 * CL_max
    
    def take_off_loading(self, CL_TO_list=None):
        """
        Calculate W/P vs W/S curves for an array of CL_TO values.
        """
        if CL_TO_list is None:
            CL_TO_list = [self.CL_max_TO / (1.1 ** 2)]
        
        W_S_range = np.linspace(1, self.x_max, 1000)
        results = []
        
        for CL_TO in CL_TO_list:
            effective_CL = CL_TO / (1.1 ** 2)
            W_P_range = (self.TOP / W_S_range) * effective_CL * (self.rho / self.rho0)
            results.append((W_P_range, W_S_range, CL_TO))
        
        return results
    
    def landing_loading(self, CL_max):
        """Calculate Wing Loading for given velocity and maximum lift coefficient in landing phase"""
        return (CL_max * self.rho * (self.landing_distance / 0.5915)) / (2 * self.f)
    
    def propellor_performance(self, prop_list=None):
        """
        Calculate W/P vs W/S curves for propellor settings.
        """
        if prop_list is None:
            prop_list = [1]
        
        W_S_range = np.linspace(1, self.x_max, 1000)
        results = []
        
        for prop in prop_list:
            W_P_range = prop * (self.n_p * (((self.rho / self.rho0) ** 0.75) * ((((self.CD0 * 0.5 * self.rho * (self.V_cruise ** 3)) / (W_S_range)) + ((W_S_range) / (np.pi * self.A * self.e * 0.5 * self.rho * self.V_cruise))) ** (-1))))
            results.append((W_P_range, W_S_range, prop))
        
        return results
    
    def plot_wing_loading_constraints(self, CL_TO_list=None):
        """Plot Wing Loading Constraints for Stall, Cruise, and Takeoff (with multiple CL values)"""
        
        plt.figure(figsize=(10, 6))
        
        # Stall constraint
        stall_lines = []
        for idx, CL in enumerate(self.CL_max_clean):
            WL_stall = self.wing_loading(self.V_stall, CL)
            stall_lines.append(WL_stall)
            color = f"C{idx % 10}"
            plt.axvline(x=WL_stall, color=color, linestyle='--', label=f"Stall CL_land={CL:.2f} @ {round(WL_stall, 2)}")
        WL_stall_limit = min(stall_lines)
        plt.axvspan(WL_stall_limit, self.x_max, color='red', alpha=0.15)
        
        # Cruise constraint
        cruise_lines = []
        for idx, CL in enumerate(self.CL_max_clean):
            WL_cruise = self.wing_loading(self.V_cruise, CL)
            cruise_lines.append(WL_cruise)
            color = f"C{idx % 10}"
            plt.axvline(x=WL_cruise, color=color, linestyle='--', label=f"Cruise CL_clean={CL:.2f} @ {round(WL_cruise, 2)}")
        WL_cruise_limit = min(cruise_lines)
        plt.axvspan(WL_cruise_limit, self.x_max, color='red', alpha=0.15)
        
        # Take-off constraint
        takeoff_curves = self.take_off_loading(CL_TO_list=self.CL_max_TO)
        takeoff_lines = []
        for idx, (W_P_range, W_S_range, CL_TO) in enumerate(takeoff_curves):
            label = f"Takeoff CL_TO={CL_TO:.2f}"
            color = f"C{idx % 10}"
            plt.plot(W_S_range, W_P_range, label=label, color=color)
            takeoff_lines.append(W_S_range[0])  # The W/S limit of the takeoff curve
        WL_takeoff_limit = min(takeoff_lines)
        plt.fill_between(W_S_range, W_P_range, self.y_max, color='red', alpha=0.15)
        
        # Landing constraint
        landing_lines = []
        for idx, CL in enumerate(self.CL_max_land):
            WL_land = self.landing_loading(CL)
            landing_lines.append(WL_land)
            color = f"C{idx % 10}"
            plt.axvline(x=WL_land, color=color, linestyle='--', label=f"Landing CL_land={CL:.2f} @ {round(WL_land, 2)}")
        WL_landing_limit = min(landing_lines)
        plt.axvspan(WL_landing_limit, self.x_max, color='red', alpha=0.15)
        
        # Propellor constraint
        propellor_curves = self.propellor_performance(prop_list=self.prop_setting)
        propellor_lines = []
        for idx, (W_P_range, W_S_range, prop) in enumerate(propellor_curves):
            label = f"Propellor Power Setting={(prop*100):.2f}%"
            color = f"C{idx % 10}"
            plt.plot(W_S_range, W_P_range, label=label, color=color)
            propellor_lines.append(W_S_range[0])  # The W/S limit of the propellor curve
        WL_propellor_limit = min(propellor_lines)
        plt.fill_between(W_S_range, W_P_range, self.y_max, color='red', alpha=0.15)
        
        # Plotting the constraints and W/P vs W/S
        plt.xlim(0, self.x_max)
        plt.ylim(0, self.y_max)
        plt.xlabel("Wing Loading (N/m²) / W/S (N/m²)")
        plt.ylabel("Power-to-Weight Ratio (W/P) [N/W]")
        plt.title("Wing Loading vs Power-to-Weight (with Multiple CL Values)")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()


# Example Usage
if __name__ == "__main__":
    aircraft = WP_WS_Diagram(
        
        # Aerodynamics
        S              = 16,    # Wing area
        A              = 7.5,   # Aspect Ratio
        e              = 0.80,  # Oswald efficiency factor
        CD0            = 0.030, # Parasite drag
        h              = 0,     # Altitude
        L              = 2,     # Chord length
        
        # Weight estimation
        W              = 6000,  # Weight
        f              = 0.9,   # Weight fraction (fuel)
        
        # Velocity
        V_cruise       = 35,    # Pre-given cruise speed
        V_stall        = 25,    # Pre-given stall speed
        
        # Design pre-set values
        TOP              = 150, # Take-off parameter (see ADSEE graph with take-off distance)
        landing_distance = 700, # Landing distance
        
        # Test variables
        CL_MAX_clean   = [1.7, 1.8, 1.9], # Possible clean CL
        CL_MAX_TO      = [1.8, 1.9, 2.0], # Possible take-off CL
        CL_MAX_land    = [1.9, 2.0, 2.1], # Possible landing CL 
        prop_setting   = [1.0, 0.9],       # Different output settings of the propellor
        
        # Propulsion
        n_p            = 0.8,  # Propellor efficiency
    )
    
    aircraft.plot_wing_loading_constraints()
