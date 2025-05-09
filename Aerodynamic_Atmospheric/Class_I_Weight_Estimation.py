import numpy as np
import matplotlib.pyplot as plt
from Requirements import *

class Class_I_Weight_Estimation():
    """
    Class I Weight Estimation based on assumed mass fractions and MTOW iteration,
    including fuel and battery usage estimation.
    """
    def __init__(self, payload_weight_kg=100,
                    residual_fuel_fraction=0.0 ,empty_weight_fraction=0.5, initial_mtow_guess_kg=1000, iteration_limit=100, tolerance=0.1, 
                    W1_WTO=1, W2_W1=1, W3_W2=1, W4_W3=1, W5_W4=1, W6_W5=1, W7_W6=1, W8_W7=1, Wfinal_W8=1,
                    n_p=0, c_p=0, g=0, A=0, e=0, C_D_0=0, V_cruise=0,
                    battery_energy_required_Wh=0, battery_specific_energy_Wh_per_kg=200):
        
        # Requirements input
        self.payload_weight_kg = payload_weight_kg
        
        # Battery inputs
        self.battery_energy_required_Wh = battery_energy_required_Wh
        self.battery_specific_energy_Wh_per_kg = battery_specific_energy_Wh_per_kg
        self.battery_mass_kg = (battery_energy_required_Wh / battery_specific_energy_Wh_per_kg) if battery_specific_energy_Wh_per_kg > 0 else 0
        
        # Inputs for first weight estimation
        self.fuel_fraction = 1 - (W1_WTO * W2_W1 * W3_W2 * W4_W3 * W5_W4 * W6_W5 * W7_W6 * W8_W7 * Wfinal_W8)
        self.residual_fuel_fraction = residual_fuel_fraction
        self.empty_weight_fraction = empty_weight_fraction
        self.MTOW_kg = initial_mtow_guess_kg
        self.iteration_limit = iteration_limit
        self.tolerance = tolerance
        
        # Fuel estimation fractions
        self.W1_WTO = W1_WTO
        self.W2_W1 = W2_W1
        self.W3_W2 = W3_W2
        self.W4_W3 = W4_W3
        self.W5_W4 = W5_W4
        self.W6_W5 = W6_W5
        self.W7_W6 = W7_W6
        self.W8_W7 = W8_W7
        self.Wfinal_W8 = Wfinal_W8
        
        # Range estimation variables
        self.n_p = n_p
        self.c_p = (c_p * 0.453592) / (745.7 * 3600)
        self.g = g
        self.A = A
        self.e = e
        self.C_D_0 = C_D_0
        self.V_cruise = V_cruise
        
        self.results = {}
        self.converged = False
        
        self.Determine_MTOW()
        self.Determine_Used_Fuel()
        self.Determine_MLW()
        #self.Determine_Brequet_Range()
        #self.Determine_Brequet_Endurance()
    
    def Determine_MTOW(self):
        for _ in range(self.iteration_limit):
            OEW = self.MTOW_kg * self.empty_weight_fraction
            fuel_weight = self.MTOW_kg * self.fuel_fraction
            new_MTOW = OEW + fuel_weight + self.payload_weight_kg + self.battery_mass_kg
            
            if abs(new_MTOW - self.MTOW_kg) < self.tolerance:
                self.converged = True
                break
            self.MTOW_kg = new_MTOW
        
        if self.converged:
            self.estimated_MTOW = self.MTOW_kg
            self.estimated_OEW = self.MTOW_kg * self.empty_weight_fraction
            self.estimated_fuel_weight = self.MTOW_kg * self.fuel_fraction
            self.results["Maximum Take-off Weight [kg]"] = round(self.estimated_MTOW, 3)
            self.results["Operation Empty Weight [kg]"] = round(self.estimated_OEW, 3)
            self.results["Fuel Weight [kg]"] = round(self.estimated_fuel_weight, 3)
            self.results["Battery Mass [kg]"] = round(self.battery_mass_kg, 3)
    
    def Determine_Used_Fuel(self):
        M_ff = (self.W1_WTO * self.W2_W1 * self.W3_W2 * self.W4_W3 *
                self.W5_W4 * self.W6_W5 * self.W7_W6 * self.W8_W7 * self.Wfinal_W8)
        self.W_f_used = (1 - M_ff) * self.MTOW_kg
        self.results["Used Fuel Estimate [kg]"] = round(self.W_f_used, 3)
    
    def Determine_MLW(self):
        if self.converged:
            residual_fuel_weight = self.estimated_fuel_weight * self.residual_fuel_fraction
            mlw = self.estimated_OEW + self.payload_weight_kg + self.battery_mass_kg + residual_fuel_weight
            self.results["Maximum Landing Weight [kg]"] = round(mlw, 3)
    
    def Determine_Maximum_Lift_Drag_Ratio(self):
        self.L_D = np.sqrt((np.pi * self.A * self.e) / (4 * self.C_D_0))
    
    def Determine_Brequet_Range(self):
        self.Determine_Maximum_Lift_Drag_Ratio()
        self.Range_Brequet = (self.n_p / (self.c_p * self.g)) * self.L_D * np.log(self.MTOW_kg / (self.MTOW_kg - self.W_f_used))
        self.results["Estimated Range [km]"] = round(self.Range_Brequet, 3)
    
    def Determine_Brequet_Endurance(self):
        self.Determine_Maximum_Lift_Drag_Ratio()
        self.Endurance_Brequet = (self.n_p / (self.V_cruise * self.c_p * self.g)) * self.L_D * np.log(self.MTOW_kg / (self.MTOW_kg - self.W_f_used))
        self.results["Estimated Endurance [s]"] = round(self.Endurance_Brequet, 3)
    
    def __str__(self):
        self.results["Estimated Weight [N]"] = round(self.estimated_MTOW * self.g, 3)
        output = ["Class I Weight Estimation Results:\n"]
        max_key_length = max(len(key) for key in self.results.keys())
        header = f"{'Parameter'.ljust(max_key_length)} | Value"
        output.append(header)
        output.append("-" * len(header))
        
        for key, value in self.results.items():
            output.append(f"{key.ljust(max_key_length)} | {value:>13,.3f}")
        
        return "\n".join(output)
    
    def Plot_Payload_Range_Diagram(self):
        MTOW = self.estimated_MTOW
        OEW = self.estimated_OEW
        Fuel = self.estimated_fuel_weight
        Battery = self.battery_mass_kg
        Payload_max = MTOW - OEW - Fuel - Battery
        MZFW = OEW + Payload_max + Battery

        range_A = 0
        payload_A = Payload_max
        range_B = self.Range_Brequet
        payload_B = Payload_max

        Range_max_fuel_only = (self.n_p / (self.c_p * self.g)) * self.L_D * np.log(MTOW / OEW)
        range_C = Range_max_fuel_only
        payload_C = 0

        plt.figure(figsize=(10,6))
        plt.plot([range_A, range_B, range_C], [payload_A, payload_B, payload_C], marker='o', label='Payload-Range Curve')
        plt.axhline(OEW, color='gray', linestyle='--', label='OEW')
        plt.axhline(MTOW, color='red', linestyle='--', label='MTOW')
        plt.axhline(self.results.get("Maximum Landing Weight [kg]", 0), color='purple', linestyle='--', label='MLW')
        plt.axhline(MZFW, color='green', linestyle='--', label='MZFW')

        plt.xlabel("Range [km]")
        plt.ylabel("Payload [kg]")
        plt.title("Payload-Range Diagram")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    req = Requirements()
    
    Class_I_Weight_Estimate = Class_I_Weight_Estimation(
        
        # Estimates and given weights
        payload_weight_kg      = req.payload_weight_kg,
        empty_weight_fraction  = 0.7,
        initial_mtow_guess_kg  = 15 * req.payload_weight_kg,
        
        # Fuel fractions
        residual_fuel_fraction = 0.00,
        W1_WTO     = 1, #0.997,
        W2_W1      = 1, #0.995,
        W3_W2      = 1, #0.996,
        W4_W3      = 1, #0.998,
        W5_W4      = 1, #0.931,
        W6_W5      = 1, #0.999,
        W7_W6      = 1, #0.998,
        W8_W7      = 1, #0.995,
        Wfinal_W8  = 1, #0.997,
        
        # Battery additions
        battery_energy_required_Wh        = 20000, # 10 kWh battery
        battery_specific_energy_Wh_per_kg = 435,   # Li-ion (200 Wh/kg)
        
        # Properties for range and endurance
        #n_p        = 0.82,
        #c_p        = 0.3,
        g          = 9.80665,
        #A          = 25,
        #e          = 0.9,
        #C_D_0      = 0.020,
        #V_cruise   = 25,
    )
    
    print(Class_I_Weight_Estimate)
    #Class_I_Weight_Estimate.Plot_Payload_Range_Diagram()
