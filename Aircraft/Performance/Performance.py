import numpy as np
import matplotlib.pyplot as plt
import aerosandbox as asb

class Performance:
    
    def __init__(self, ac, plot):
        self.ac = ac
        self.show_plot = plot
        
        # Inputs from aircraft inputs and weight estimation
        self.n_p = ac.inputs.n_p
        self.altitude = ac.inputs.h_cruise
        self.isa = ac.ISA
        self.rho = ac.rho
        self.g = ac.inputs.g
        self.estimated_MTOM = ac.class_I.estimated_MTOM
        self.weight_N = self.estimated_MTOM * self.g
        self.battery_power_available = ac.inputs.battery_power_available  # in W
        self.battery_specific_energy_Wh_per_kg = ac.inputs.battery_specific_energy_Wh_per_kg
        
        # Propulsion choice and duration
        self.propulsion_type = ac.inputs.propulsion_type.lower()  # 'battery' or 'fuel'
        self.desired_endurance_sec = ac.inputs.endurance * 24 * 3600  # days → seconds
        
        # Run functions
        self.compute_drag_and_power()
        if self.propulsion_type == 'battery':
            self.handle_battery_endurance()
        elif self.propulsion_type == 'hydrogen':
            self.handle_hydrogen_endurance()
        else:
            raise ValueError("Unknown propulsion type. Use 'battery' or 'fuel'.")
    
    def compute_drag_and_power(self):
    
        airplane = self.ac.configuration
        S = airplane.s_ref
        W = self.ac.class_I.estimated_MTOW
        
        velocities = np.linspace(1, 20, 380)  # avoid V=0
        V_list = []
        DV_list = []
        can_fly_mask = []
        
        for V in velocities:
            # Set up properties
            AoAs = np.linspace(-2, 15, 85)
            op_point = asb.OperatingPoint(velocity=V, alpha=AoAs, beta=0)
            aero = asb.AeroBuildup(airplane=airplane, op_point=op_point).run()
            
            # Get properties from aerodynamics and guess the trim angle
            CL_array = aero['CL']
            L_array = 0.5 * self.rho * V**2 * S * CL_array
            diff_array = np.abs(L_array - W)
            best_idx = np.argmin(diff_array)
            trim_guess = AoAs[best_idx]
            
            # Find the exact trim angle
            AoAs_fine = np.linspace(trim_guess - 1.0, trim_guess + 1.0, 100)
            op_point_fine = asb.OperatingPoint(velocity=V, alpha=AoAs_fine, beta=0)
            aero_fine = asb.AeroBuildup(airplane=airplane, op_point=op_point_fine).run()
            
            # Define the trim condition
            CL_fine = aero_fine['CL']
            CD_fine = aero_fine['CD']
            L_fine = 0.5 * self.rho * V**2 * S * CL_fine
            diff_fine = np.abs(L_fine - W)
            valid_lift_mask = L_fine >= W # Find the first valid AoA where L >= W 
            
            # Check if the L >= W
            if np.any(valid_lift_mask):
                best_idx_fine = np.argmax(valid_lift_mask) 
                best_CD = CD_fine[best_idx_fine]
                best_L = L_fine[best_idx_fine]
            else:
                best_idx_fine = np.argmin(diff_fine)
                best_CD = CD_fine[best_idx_fine]
                best_L = L_fine[best_idx_fine]
            
            # Find drag and power for trim condition
            D = 0.5 * self.rho * V**2 * S * best_CD
            DV = D * V
            V_list.append(V)
            DV_list.append(DV)
            if best_L >= W:
                print(f"Trim succeeded at V = {V:.2f} m/s and found Power Required of P_req = {(D*V):.2f} W")
                can_fly_mask.append(True)
            else:
                print(f"Trim failed at V = {V:.2f} m/s")
                can_fly_mask.append(False)
        
        # Collect calculation results
        V_array = np.array(V_list)
        DV_array = np.array(DV_list)
        can_fly_mask = np.array(can_fly_mask)
        if not np.any(can_fly_mask):
            raise RuntimeError("No velocity found where lift exceeds weight.")
        
        # Define valid optimum space (no stall i.e.)
        V_valid = V_array[can_fly_mask]
        DV_valid = DV_array[can_fly_mask]
        
        # Find stall and optimum for minimum power required
        self.stall_speed = V_valid[0]  # First flyable speed
        min_power_idx = np.argmin(DV_valid)
        self.optimum_V = V_valid[min_power_idx]
        self.battery_power_required = DV_valid[min_power_idx]
        
        plt.figure(figsize=(10, 6))
        plt.plot(V_array[~can_fly_mask], DV_array[~can_fly_mask], 'x', color='gray', label='Below Stall (Not Flyable)')
        plt.plot(V_array[can_fly_mask], DV_array[can_fly_mask], '-', color='blue', label='Flyable Power')
        plt.axvline(self.stall_speed, color='red', linestyle='--', label='Stall Speed')
        plt.axvline(self.optimum_V, color='green', linestyle='--', label='Min Power Speed')
        plt.xlabel("Velocity [m/s]")
        plt.ylabel("Power Required [W]")
        plt.title("Power Required vs. Velocity")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        if self.show_plot:
            plt.show()
    
    def handle_battery_endurance(self):
        """Battery-specific endurance and energy calculation."""
        energy_required_Wh = self.battery_power_required * self.desired_endurance_sec / 3600
        battery_mass_kg = energy_required_Wh / self.battery_specific_energy_Wh_per_kg
        self.battery_mass_kg = battery_mass_kg
    
    def handle_hydrogen_endurance(self):
        """Estimate required hydrogen mass for given endurance and power."""
        # Find required energy
        hydrogen_specific_energy_Wh_per_kg = self.ac.inputs.hydrogen_specific_energy_Wh_per_kg  
        energy_required_Wh = self.battery_power_required * self.desired_endurance_sec / 3600 
        
        # Find required mass of hydrogen
        hydrogen_mass_kg = energy_required_Wh / hydrogen_specific_energy_Wh_per_kg
        self.hydrogen_mass_kg = hydrogen_mass_kg
        
        # Find required tank size
        hydrogen_density_kg_per_m3 = self.ac.inputs.hydrogen_density_kg_per_m3 
        self.hydrogen_tank_volume_m3 = self.hydrogen_mass_kg / hydrogen_density_kg_per_m3
        
        # Find required tank mass
        tank_mass_fraction = self.ac.inputs.tank_mass_fraction  
        self.hydrogen_tank_mass_kg = self.hydrogen_mass_kg * tank_mass_fraction
    
    def __str__(self):
        output = ["\n"]
        output.append("Performance Estimation Results:")
        results_data = {
            "Optimum Velocity for Endurance [m/s]": getattr(self, "optimum_V", None),
            "Found Stall Velocity while in trim [m/s]": getattr(self, "stall_speed", None),
            "Battery Power Required for required endurance [W]": getattr(self, "battery_power_required", None),
            "Battery Mass Required for required endurance [kg]": getattr(self, "battery_mass_kg", None),
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
