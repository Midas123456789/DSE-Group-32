import numpy as np
import matplotlib.pyplot as plt
from Requirements import *

class Class_I_Weight_Estimation():
    """
    Class I Weight Estimation based on assumed mass fractions and MTOW iteration,
    including fuel and battery usage estimation.
    """
    def __init__(self, payload_weight_kg,
                    residual_fuel_fraction ,empty_weight_fraction, initial_mtow_guess_kg, iteration_limit, tolerance, 
                    W1_WTO, W2_W1, W3_W2, W4_W3, W5_W4, W6_W5, W7_W6, W8_W7, Wfinal_W8,
                    n_p, c_p, g, A, e, CD0,
                    battery_power_available, battery_specific_energy_Wh_per_kg, endurance):
        
        # Requirements input
        self.payload_weight_kg = payload_weight_kg
        
        # Battery inputs
        self.battery_power_available = battery_power_available
        self.battery_specific_energy_Wh_per_kg = battery_specific_energy_Wh_per_kg
        self.battery_mass_kg = ((battery_power_available * (endurance * 24)) / battery_specific_energy_Wh_per_kg) if battery_specific_energy_Wh_per_kg > 0 else 0
        
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
        self.CD0 = CD0
        
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
            self.estimated_MTOM = self.MTOW_kg
            self.estimated_OEM = self.MTOW_kg * self.empty_weight_fraction
            self.estimated_fuel_mass = self.MTOW_kg * self.fuel_fraction
            
            self.estimated_MTOW = self.MTOW_kg * self.g
            self.estimated_OEW = self.MTOW_kg * self.empty_weight_fraction * self.g
            self.estimated_fuel_weight = self.MTOW_kg * self.fuel_fraction * self.g
            
            self.results["Maximum Take-off Weight [kg]"] = round(self.estimated_MTOM, 3)
            self.results["Estimated Weight [N]"] = round(self.estimated_MTOM * self.g, 3)
            self.results["Operating Empty Weight [kg]"] = round(self.estimated_OEM, 3)
            self.results["Fuel Weight [kg]"] = round(self.estimated_fuel_mass, 3)
            self.results["Battery Mass [kg]"] = round(self.battery_mass_kg, 3)
        else:
            print(50*"-")
            print("Stopped in Class I, MTOW diverged!!")
            print(50*"-") 
    
    def Determine_Used_Fuel(self):
        M_ff = (self.W1_WTO * self.W2_W1 * self.W3_W2 * self.W4_W3 *
                self.W5_W4 * self.W6_W5 * self.W7_W6 * self.W8_W7 * self.Wfinal_W8)
        self.W_f_used = (1 - M_ff) * self.MTOW_kg
        self.results["Used Fuel Estimate [kg]"] = round(self.W_f_used, 3)
    
    def Determine_MLW(self):
        if self.converged:
            residual_fuel_weight = self.estimated_fuel_mass * self.residual_fuel_fraction
            mlw = self.estimated_OEW + self.payload_weight_kg + self.battery_mass_kg + residual_fuel_weight
            self.results["Maximum Landing Weight [kg]"] = round(mlw, 3)
    
    def Determine_Maximum_Lift_Drag_Ratio(self):
        self.L_D = np.sqrt((np.pi * self.A * self.e) / (4 * self.CD0))
    
    def __str__(self):
        output = ["\n"]
        output.append("Class I Weight Estimation Results:")
        
        # Prepare list of key-value pairs to print
        results_data = {
            "Maximum Take-off Mass [kg]": getattr(self, "estimated_MTOM", None),
            "Maximum Take-off Weight [N]": getattr(self, "estimated_MTOW", None),
            "Operating Empty Mass [kg]": getattr(self, "estimated_OEM", None),
            "Operating Empty Weight [N]": getattr(self, "estimated_OEW", None),
            "Fuel Mass [kg]": getattr(self, "estimated_fuel_mass", None),
            "Fuel Weight [N]": getattr(self, "estimated_fuel_weight", None),
            "Battery Mass [kg]": self.battery_mass_kg,
            "Used Fuel Estimate [kg]": getattr(self, "W_f_used", None),
            "Maximum Landing Weight [kg]": self.results.get("Maximum Landing Weight [kg]", None)
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
        battery_power_available           = 20000, # 10 kWh battery
        battery_specific_energy_Wh_per_kg = 435,   # Li-ion (200 Wh/kg)
        
        # Properties for range and endurance
        n_p        = 0.82,
        c_p        = 0.3,
        g          = 9.80665,
        A          = 25,
        e          = 0.9,
        CD0        = 0.020,
    )
    
    print(Class_I_Weight_Estimate)
    #Class_I_Weight_Estimate.Plot_Payload_Range_Diagram()