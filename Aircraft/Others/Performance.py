import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

class Performance:
    def __init__(self, ac):
        self.ac = ac
        # Inputs
        self.A = ac.inputs.A
        self.e = ac.inputs.e
        self.CD0 = ac.inputs.CD0
        self.n_p = ac.inputs.n_p
        self.S = ac.inputs.S
        self.altitude = ac.inputs.h_cruise
        self.isa = ac.ISA
        self.rho = self.isa.results[self.altitude]["Density [kg/m3]"]
        self.g = self.isa.results[self.altitude]["Gravity [m/s2]"]
        self.estimated_MTOM = ac.class_I.estimated_MTOM
        self.weight_N = self.estimated_MTOM * self.g
        self.L_D_endurance = np.sqrt((3 * np.pi * self.A * self.e) / self.CD0)  # Maximum L/D ratio for endurance
        self.desired_endurance = ac.inputs.endurance * 24 * 60 * 60  # days to seconds
        self.battery_power_available = ac.inputs.battery_power_available  # kW to W
        self.design_for_endurance()

    def compute_drag_and_power(self, V):
        # Compute drag forces and power required for a given velocity V
        q = 0.5 * self.rho * V**2
        parasite_drag = q * self.S * self.CD0
        induced_drag = (self.weight_N**2) / (q * self.S * np.pi * self.A * self.e)
        total_drag = parasite_drag + induced_drag
        power_required = total_drag * V / self.n_p  # Power in Watts
        return total_drag, power_required
    
    def design_for_endurance(self):
        
        # Minimize power required and find optimal velocity
        result = minimize_scalar(lambda V: self.compute_drag_and_power(V)[1], bounds=(10, 200), method='bounded')
        optimal_velocity = result.x
        power_required_minimum = result.fun
        
        # Calculate drag forces at optimal velocity
        total_drag, power_required_minimum = self.compute_drag_and_power(optimal_velocity)
        energy_required_Wh = power_required_minimum * (self.desired_endurance / 3600)
        battery_mass_kg = energy_required_Wh / self.ac.inputs.battery_specific_energy_Wh_per_kg
        
        self.optimum_V = optimal_velocity
        self.battery_power_required = power_required_minimum
        self.battery_mass_kg = battery_mass_kg
    
    def plot_drag_to_power_for_max_endurance(self):
        # Visualize drag and power required
        velocities = np.linspace(0, 200, 200)
        drag_values, power_values = [], []
        
        for V in velocities:
            total_drag, power_required = self.compute_drag_and_power(V)
            drag_values.append(total_drag)
            power_values.append(power_required)
        
        # Plot drag and power required separately for visualization
        plt.figure(figsize=(8, 6))
        plt.plot(velocities, drag_values, label='Total Drag (D)', color='b')
        plt.plot(velocities, power_values, label='Power Required (D * V)', color='purple')
        plt.title('Velocity vs Drag & Power Required')
        plt.xlabel('Velocity (m/s)')
        plt.ylabel('Drag (N) & Power (W)')
        plt.grid(True)
        plt.legend()
        plt.show()