import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

class Performance:
    def __init__(self, ac):
        self.ac = ac

        # Inputs from aircraft
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
        self.L_D_endurance = np.sqrt((3 * np.pi * self.A * self.e) / self.CD0)

        # Propulsion choice and duration
        self.propulsion_type = ac.inputs.propulsion_type.lower()  # 'battery' or 'fuel'
        self.desired_endurance_sec = ac.inputs.endurance * 24 * 3600  # days → seconds

        self.battery_power_available = ac.inputs.battery_power_available  # in W
        self.battery_specific_energy_Wh_per_kg = ac.inputs.battery_specific_energy_Wh_per_kg

        self.design_for_endurance()

    def design_for_endurance(self):
        """Optimizes flight velocity for minimum power (max endurance) and computes energy/fuel needs."""
        # Find velocity that minimizes power required
        result = minimize_scalar(
            lambda V: self.ac.aero.compute_drag_and_power(V)[1],
            bounds=(10, 200),
            method='bounded'
        )
        self.optimum_V = result.x
        _, self.battery_power_required = self.ac.aero.compute_drag_and_power(self.optimum_V)

        if self.propulsion_type == 'battery':
            self.handle_battery_endurance()
        elif self.propulsion_type == 'hydrogen':
            self.handle_hydrogen_endurance()
        else:
            raise ValueError("Unknown propulsion type. Use 'battery' or 'fuel'.")

    def handle_battery_endurance(self):
        """Battery-specific endurance and energy calculation."""
        energy_required_Wh = self.battery_power_required * self.desired_endurance_sec / 3600
        battery_mass_kg = energy_required_Wh / self.battery_specific_energy_Wh_per_kg

        self.battery_mass_kg = battery_mass_kg

    def handle_hydrogen_endurance(self):
        """Estimate required hydrogen mass for given endurance and power."""
        hydrogen_specific_energy_Wh_per_kg = 33333  # 120 MJ/kg in Wh
        energy_required_Wh = self.battery_power_required * self.desired_endurance_sec / 3600  # in Wh

        hydrogen_mass_kg = energy_required_Wh / hydrogen_specific_energy_Wh_per_kg

        self.hydrogen_mass_kg = hydrogen_mass_kg
        
        self.estimate_hydrogen_tank_size_and_weight()
    
    def estimate_hydrogen_tank_size_and_weight(self):
        """
        Estimate hydrogen tank volume and mass based on storage assumptions.
        Assumes compressed hydrogen gas at 700 bar and 20°C.
        """
        # Constants for compressed hydrogen at 700 bar and 288 K
        hydrogen_density_kg_per_m3 = 42  # Approximate density at 700 bar, 15°C
        
        # Volume = mass / density
        tank_volume_m3 = self.hydrogen_mass_kg / hydrogen_density_kg_per_m3
        
        # Estimate tank mass based on mass fraction or specific factor
        tank_mass_fraction = 0.15  # e.g., 15% of hydrogen mass
        tank_mass_kg = self.hydrogen_mass_kg * tank_mass_fraction

        # Store in object
        self.hydrogen_tank_volume_m3 = tank_volume_m3
        self.hydrogen_tank_mass_kg = tank_mass_kg


    def plot_drag_to_power_for_max_endurance(self):
        """Visualizes drag and power required across velocities."""
        velocities = np.linspace(5, 100, 100)
        drag_values, power_values = [], []

        for V in velocities:
            total_drag, power_required = self.compute_drag_and_power(V)
            drag_values.append(total_drag)
            power_values.append(power_required)

        plt.figure(figsize=(8, 6))
        plt.plot(velocities, drag_values, label='Total Drag (N)', color='b')
        plt.plot(velocities, power_values, label='Power Required (W)', color='purple')
        plt.title('Velocity vs Drag & Power Required')
        plt.xlabel('Velocity (m/s)')
        plt.ylabel('Drag (N) / Power (W)')
        plt.grid(True)
        plt.legend()
        plt.show()

    def __str__(self):
        output = ["\n"]
        output.append("Performance Estimation Results:")
        results_data = {
            "Optimum Velocity for Endurance [m/s]": getattr(self, "optimum_V", None),
            "Battery Power Required [W]": getattr(self, "battery_power_required", None),
            "Battery Mass Required [kg]": getattr(self, "battery_mass_kg", None),
            "Hydrogen Mass Required [kg]": getattr(self, "hydrogen_mass_kg", None),
            "Hydrogen Tank Volume [m³]": getattr(self, "hydrogen_tank_volume_m3", None),
            "Hydrogen Tank Mass [kg]": getattr(self, "hydrogen_tank_mass_kg", None),
        }

        max_key_length = max(len(k) for k in results_data)
        header = f"{'Parameter'.ljust(max_key_length)} | Value"
        output.append(header)
        output.append("-" * len(header))

        for key, value in results_data.items():
            if value is not None:
                output.append(f"{key.ljust(max_key_length)} | {value:>13,.3f}")
            else:
                output.append(f"{key.ljust(max_key_length)} | {'Not computed':>13}")

        return "\n".join(output)
