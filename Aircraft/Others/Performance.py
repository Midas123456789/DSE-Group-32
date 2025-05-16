import numpy as np
from ISA_Calculator import *
import matplotlib.pyplot as plt


class Performance:
    def __init__(self, ac, plot=False):
        self.ac = ac

        # Inputs
        self.A = ac.inputs.A
        self.e = ac.inputs.e
        self.CD0 = ac.inputs.CD0
        self.n_p = ac.inputs.n_p
        self.c_p = (ac.inputs.c_p * 0.453592) / (745.7 * 3600)
        self.S = ac.inputs.S
        self.altitude = ac.inputs.h_cruise

        # ISA
        self.isa = ac.ISA
        self.rho = self.isa.results[self.altitude]["Density [kg/m3]"]
        self.g = self.isa.results[self.altitude]["Gravity [m/s2]"]

        # Class I estimates
        self.estimated_MTOM = ac.class_I.estimated_MTOM

        self.W_f_used = ac.class_I.W_f_used  # Assuming all fuel used for endurance
        self.battery_power_available = ac.inputs.battery_power_available

        # Constants
        self.weight_N = self.estimated_MTOM * self.g
        self.L_D_endurance = np.sqrt((3 * np.pi * self.A * self.e) / self.CD0) # Maximum endurance
        
        self.plot = plot
        self.plot_endurance_vs_velocity()

    def plot_endurance_vs_velocity(self, plot=False):
        V_range = np.linspace(10, 150, 1500)

        total_endurance_list = []
        fuel_endurance_list = []
        battery_endurance_list = []
        
        for V in V_range:
            # Fuel-based endurance
            if self.W_f_used > 0:
                fuel_endurance = (self.n_p / (V * self.c_p * self.g)) * self.L_D_endurance * np.log(self.estimated_MTOM / (self.estimated_MTOM - self.W_f_used))
            else:
                fuel_endurance = 0

            # Battery-based endurance
            if self.battery_power_available > 0:
                q = 0.5 * self.rho * V**2
                parasite_drag = q * self.S * self.CD0
                induced_drag = (self.weight_N**2) / (q * self.S * np.pi * self.A * self.e)
                total_drag = parasite_drag + induced_drag
                power_required_W = total_drag * V / self.n_p
                battery_endurance = ((self.battery_power_available * 3600) / power_required_W)  # [seconds]
            else:
                battery_endurance = 0

            total_endurance = fuel_endurance + battery_endurance

            total_endurance_list.append(total_endurance)
            fuel_endurance_list.append(fuel_endurance)
            battery_endurance_list.append(battery_endurance)

        # Plotting
        plt.figure(figsize=(10, 6))
        plt.plot(V_range, total_endurance_list, label='Total Endurance', color='black')
        plt.plot(V_range, fuel_endurance_list, label='Fuel Endurance', linestyle='--', color='blue')
        plt.plot(V_range, battery_endurance_list, label='Battery Endurance', linestyle='--', color='green')

        plt.xlabel("Velocity [m/s]")
        plt.ylabel("Endurance [s]")
        plt.title(f"Endurance vs Velocity (Battery {self.battery_power_available}[Wh], Fuel {self.W_f_used}[kg])")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()

        max_idx = np.argmax(total_endurance_list)
        self.optimum_V = V_range[max_idx]
        self.max_endurance = max(total_endurance_list)

        self.plot = plot
        if self.plot:
            plt.show()


if __name__ == "__main__":
    # Call the function with the necessary inputs
    per = Performance(ac)
    
    ac = Aircraft(
        S = 30,
        A = 25,
        e = 0.8,
        CD0 = 0.040,
        n_p = 0.95,
        c_p = (0.6 * 0.453592) / (745.7 * 3600), # [kg/Ws]
        g = 9.80665,
        estimated_MTOM = 500,
        W_f_used = 0,
        battery_power_available = 20000, # [Wh]
        altitude=0
    )
