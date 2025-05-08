from ISA_Calculator import *
import numpy as np
import matplotlib.pyplot as plt

class WP_WS_Diagram:
    
    def __init__(self, V_cruise=0, V_stall=0, A=0, e=0, CD0=0, h=0,
            CL_MAX_clean=None, CL_MAX_land=0, CL_MAX_TO=None, prop_setting=0, TOP=0, landing_distance=0, f=1, n_p=1):
        # Ensure clean and TO CLs are lists
        self.x_max = 1000
        self.y_max = 1
        
        self.CL_max_clean = CL_MAX_clean if isinstance(CL_MAX_clean, list) else [CL_MAX_clean]
        self.CL_max_land  = CL_MAX_land if isinstance(CL_MAX_land, list) else [CL_MAX_land]
        self.CL_max_TO    = CL_MAX_TO if isinstance(CL_MAX_TO, list) else [CL_MAX_TO]
        self.prop_setting = prop_setting if isinstance(prop_setting, list) else [prop_setting]
        
        self.V_cruise = V_cruise
        self.V_stall = V_stall
        
        self.A = A 
        self.e = e 
        self.CD0 = CD0
        
        self.f = f
        
        self.TOP = TOP  # Takeoff Parameter (depends on aircraft class)
        self.landing_distance = landing_distance
        self.n_p = n_p
        
        self.h = h if isinstance(h, list) else [h]
        self.rho = []
        self.rho0 = 1.225
        self.set_densities()
    
    def set_densities(self):
        for h in self.h:
            isa = ISA_Calculator(self.h)
            self.rho.append(isa.results[h]["Density [kg/m³]"])
    
    def wing_loading(self, V, CL_max, rho=1.225):
        """Calculate Wing Loading for given velocity and maximum lift coefficient in cruise"""
        return 0.5 * rho* V**2 * CL_max
    
    def take_off_loading(self, CL_TO_list=None, rho=1.225):
        """
        Calculate W/P vs W/S curves for an array of CL_TO values.
        """
        if CL_TO_list is None:
            CL_TO_list = [self.CL_max_TO / (1.1 ** 2)]
        
        W_S_range = np.linspace(1, self.x_max, 1000)
        results = []
        
        for CL_TO in CL_TO_list:
            effective_CL = CL_TO / (1.1 ** 2)
            W_P_range = (self.TOP / W_S_range) * effective_CL * (rho/ self.rho0)
            results.append((W_P_range, W_S_range, CL_TO))
        
        return results
    
    def landing_loading(self, CL_max, rho=1.225):
        """Calculate Wing Loading for given velocity and maximum lift coefficient in landing phase"""
        return (CL_max * rho* (self.landing_distance / 0.5915)) / (2 * self.f)
    
    def propellor_performance(self, prop_list=None, rho=1.225):
        """
        Calculate W/P vs W/S curves for propellor settings.
        """
        if prop_list is None:
            prop_list = [1]
        
        W_S_range = np.linspace(1, self.x_max, 1000)
        results = []
        
        for prop in prop_list:
            W_P_range = prop * (self.n_p * (((rho/ self.rho0) ** 0.75) * ((((self.CD0 * 0.5 * rho* (self.V_cruise ** 3)) / (W_S_range)) + ((W_S_range) / (np.pi * self.A * self.e * 0.5 * rho* self.V_cruise))) ** (-1))))
            results.append((W_P_range, W_S_range, prop))
        
        return results
    
    def plot_wing_loading_constraints(self, CL_TO_list=None):
        """Plot Wing Loading Constraints for Stall, Cruise, Takeoff, etc. across up to 4 altitudes as subplots."""

        num_altitudes = min(len(self.rho), 3)
        fig, axs = plt.subplots(1, num_altitudes, figsize=(6 * num_altitudes, 6), sharey=True)
        axs = np.atleast_1d(axs)  # Ensure axs is iterable even for 1 subplot

        for idh, rho in enumerate(self.rho[:num_altitudes]):  # Limit to 4 plots
            ax = axs[idh]

            # Stall constraint
            stall_lines = []
            for idx, CL in enumerate(self.CL_max_clean):
                WL_stall = self.wing_loading(self.V_stall, CL, rho=rho)
                stall_lines.append(WL_stall)
                color = f"C{idx % 10}"
                ax.axvline(x=WL_stall, color=color, linestyle='--', label=f"Stall CL={CL:.2f} @ {round(WL_stall, 2)}")
            WL_stall_limit = min(stall_lines)
            ax.axvspan(WL_stall_limit, self.x_max, color='red', alpha=0.15)

            # Cruise constraint
            cruise_lines = []
            for idx, CL in enumerate(self.CL_max_clean):
                WL_cruise = self.wing_loading(self.V_cruise, CL, rho=rho)
                cruise_lines.append(WL_cruise)
                color = f"C{idx % 10 + 3}"
                ax.axvline(x=WL_cruise, color=color, linestyle='--', label=f"Cruise CL={CL:.2f} @ {round(WL_cruise, 2)}")
            WL_cruise_limit = min(cruise_lines)
            ax.axvspan(WL_cruise_limit, self.x_max, color='red', alpha=0.15)

            # Take-off constraint
            takeoff_curves = self.take_off_loading(CL_TO_list=self.CL_max_TO, rho=rho)
            for idx, (W_P_range, W_S_range, CL_TO) in enumerate(takeoff_curves):
                label = f"Takeoff CL_TO={CL_TO:.2f}"
                color = f"C{idx % 10}"
                ax.plot(W_S_range, W_P_range, label=label, color=color)
            ax.fill_between(W_S_range, W_P_range, self.y_max, color='red', alpha=0.15)

            # Landing constraint
            for idx, CL in enumerate(self.CL_max_land):
                WL_land = self.landing_loading(CL, rho=rho)
                color = f"C{idx % 10 + 6}"
                ax.axvline(x=WL_land, color=color, linestyle='--', label=f"Landing CL={CL:.2f} @ {round(WL_land, 2)}")
            ax.axvspan(WL_land, self.x_max, color='red', alpha=0.15)

            # Propellor constraint
            propellor_curves = self.propellor_performance(prop_list=self.prop_setting, rho=rho)
            for idx, (W_P_range, W_S_range, prop) in enumerate(propellor_curves):
                label = f"Prop. Setting={prop*100:.0f}%"
                color = f"C{idx % 10}"
                ax.plot(W_S_range, W_P_range, label=label, color=color)
            ax.fill_between(W_S_range, W_P_range, self.y_max, color='red', alpha=0.15)

            # Axis setup
            ax.set_xlim(0, self.x_max)
            ax.set_ylim(0, self.y_max)
            ax.set_xlabel("Wing Loading W/S (N/m²)")
            if idh == 0:
                ax.set_ylabel("Weight-to-Power W/P (N/W)")
            ax.set_title(f"Altitude: {self.h[idh]} m")
            ax.grid(True)
            ax.legend(fontsize='small', loc='upper right')

        plt.suptitle("Wing Loading vs Weight-to-Power Across Altitudes", fontsize=16)
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        plt.show()



# Example Usage
if __name__ == "__main__":
    aircraft = WP_WS_Diagram(
        
        # Aerodynamics
        A              = 12,   # Aspect Ratio
        e              = 0.80,  # Oswald efficiency factor
        CD0            = 0.030, # Parasite drag
        h              = [0,10000,20000],     # Altitude
        
        # Weight estimation
        f              = 1,   # Weight fraction (fuel)
        
        # Velocity
        V_cruise       = 15,    # Pre-given cruise speed
        V_stall        = 12,    # Pre-given stall speed
        
        # Design pre-set values
        TOP              = 150, # Take-off parameter (see ADSEE graph with take-off distance)
        landing_distance = 700, # Landing distance
        
        # Test variables
        CL_MAX_clean   = [1.7], # Possible clean CL
        CL_MAX_TO      = [1.9], # Possible take-off CL
        CL_MAX_land    = [2.2], # Possible landing CL 
        prop_setting   = [1.0], # Different output settings of the propellor
        
        # Propulsion
        n_p            = 0.95,  # Propellor efficiency
    )
    
    aircraft.plot_wing_loading_constraints()
