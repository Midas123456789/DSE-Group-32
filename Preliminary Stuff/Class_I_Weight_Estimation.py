import numpy as np
from Requirements import *

class Class_I_Weight_Estimation():
    """
    Class I Weight Estimation based on assumed mass fractions and MTOW iteration,
    including fuel usage estimation.
    """
    def __init__(self, payload_weight_kg=100,
                    fuel_fraction=0.0, empty_weight_fraction=0.5, initial_mtow_guess_kg=1000, iteration_limit=100, tolerance=1.0, 
                    W1_WTO=1, W2_W1=1, W3_W2=1, W4_W3=1, W5_W4=1, W6_W5=1, W7_W6=1, W8_W7=1, Wfinal_W8=1,
                    n_p=1, c_p=100, g=10, A=1, e=1, C_D_0=0, V_cruise=10):
        
        # Requirements input
        self.payload_weight_kg      = payload_weight_kg  # Weight of payload
        
        # Inputs for first weight estimation
        self.fuel_fraction          = fuel_fraction          # Fraction of fuel compared to MTOW
        self.empty_weight_fraction  = empty_weight_fraction  # Fraction of empty weight compared to MTOW
        self.MTOW_kg                = initial_mtow_guess_kg  # Initial guessed MTOW
        self.iteration_limit        = iteration_limit        # Maximum weight estimation iterations
        self.tolerance              = tolerance              # Tolerance for convergence calling
        
        # Fuel estimation fractions
        self.W1_WTO     = W1_WTO     # Pre-Flight Operations
        self.W2_W1      = W2_W1      # Taxiing
        self.W3_W2      = W3_W2      # Take-off
        self.W4_W3      = W4_W3      # Climb
        self.W5_W4      = W5_W4      # Cruise
        self.W6_W5      = W6_W5      # Descent
        self.W7_W6      = W7_W6      # Landing
        self.W8_W7      = W8_W7      # Taxiing
        self.Wfinal_W8  = Wfinal_W8  # Post-Flight Operations
        
        # Range estimation variables
        self.n_p      = n_p       # Propeller efficiency [-]
        self.c_p      = c_p       # Specific fuel consumption [kg/J]
        self.g        = g         # Acceleration due to gravity [m/s^2]
        self.A        = A         # Wing aspect ratio [-]
        self.e        = e         # Oswald efficieny factor [-]
        self.C_D_0    = C_D_0     # Parasite drag coefficient [-]
        self.V_cruise = V_cruise  # Cruise velocity [m/s]
        
        # Dictionary to store results
        self.results = {}
        
        self.Determine_MTOW()
        self.Determine_Used_Fuel()
        self.Determine_Brequet_Range()
        self.Determine_Brequet_Endurance()
    
    def Determine_MTOW(self):
        """
        Iteratively estimates the MTOW, OEW, and fuel weight until convergence.
        """
        for _ in range(self.iteration_limit):
            OEW = self.MTOW_kg * self.empty_weight_fraction
            fuel_weight = self.MTOW_kg * self.fuel_fraction
            new_MTOW = OEW + fuel_weight + self.payload_weight_kg
            
            if abs(new_MTOW - self.MTOW_kg) < self.tolerance:
                self.converged = True
                break
            self.MTOW_kg = new_MTOW
        
        if self.converged:
            self.results["Maximum Take-off Weight [kg]"] = round(self.MTOW_kg, 3)
            self.results["Operation Empty Weight [kg]"] = round(self.MTOW_kg * self.empty_weight_fraction, 3)
            self.results["Fuel Weight [kg]"] = round(self.MTOW_kg * self.fuel_fraction, 3)
    
    def Determine_Used_Fuel(self):
        """
        Estimates the fuel used based on weight fractions during the flight phases.
        """
        M_ff = (self.W1_WTO * self.W2_W1 * self.W3_W2 * self.W4_W3 * self.W5_W4 * self.W6_W5 * self.W7_W6 * self.W8_W7 * self.Wfinal_W8)
        self.W_f_used = (1 - M_ff) * self.MTOW_kg
        self.results["Used Fuel Estimate [kg]"] = round(self.W_f_used, 3)
    
    def Determine_Maximum_Lift_Drag_Ratio(self):
        """
        Calculates the L/D ratio based on maximum approximation of a propellor aircraft.
        """
        self.L_D = np.sqrt((np.pi * self.A * self.e) / (4 * self.C_D_0))
    
    def Determine_Brequet_Range(self):
        """ 
        Estimates the Range based on fuel, weight and aerodynamic properties.
        """
        self.Determine_Maximum_Lift_Drag_Ratio()
        self.Range_Brequet = (self.n_p / (self.c_p * self.g)) * (self.L_D) * np.log(self.MTOW_kg / (self.MTOW_kg - self.W_f_used))
        self.results["Estimated Range [km]"] = round(self.Range_Brequet, 3)
    
    def Determine_Brequet_Endurance(self):
        """ 
        Estimates the Endurance based on fuel, weight and aerodynamic properties.
        """
        self.Determine_Maximum_Lift_Drag_Ratio()
        self.Endurance_Brequet = (self.n_p / (self.V_cruise * self.c_p * self.g)) * (self.L_D) * np.log(self.MTOW_kg / (self.MTOW_kg - self.W_f_used))
        self.results["Estimated Endurance [s]"] = round(self.Range_Brequet, 3)
    
    def __str__(self):
        output = ["Class I Weight Estimation Results:"]
        for key, value in self.results.items():
            output.append(f"  {key}: {value}")
        return "\n".join(output)


if __name__ == "__main__":
    req = Requirements()
    
    # Create a Class_I_Weight_Estimation object with the required parameters
    Class_I_Weight_Estimate = Class_I_Weight_Estimation(
        
        # Weight Estimation
        fuel_fraction          = 0.05,                       # Fraction of fuel compared to MTOW
        empty_weight_fraction  = 0.8,                        # Fraction of empty weight compared to MTOW
        initial_mtow_guess_kg  = 8 * req.payload_weight_kg,  # Initial guessed MTOW
        iteration_limit        = 100,                        # Maximum weight estimation iterations
        tolerance              = 1,                          # Tolerance for convergence calling
        
        # Fuel Estimation **(Data taken from Image at agriculture)**
        W1_WTO                 = 0.997,  # Pre-Flight Operations
        W2_W1                  = 0.995,  # Taxiing
        W3_W2                  = 0.996,  # Take-off
        W4_W3                  = 0.998,  # Climb
        W5_W4                  = 0.931,  # Cruise
        W6_W5                  = 0.999,  # Descent
        W7_W6                  = 0.998,  # Landing
        W8_W7                  = 0.995,  # Taxiing
        Wfinal_W8              = 0.997,  # Post-Flight Operations
        
        # Range Estimation **(Data taken from Image at agriculture)**
        n_p                    = 0.82,    # Propeller efficiency [-]
        c_p                    = 0.6,     # Specific fuel consumption [lbs/hr/hp]
        g                      = 9.80665, # Acceleration due to gravity [m/s^2]
        A                      = 12,      # Wing aspect ratio [-]
        e                      = 0.7,     # Oswald efficieny factor [-]
        C_D_0                  = 0.060,   # Parasite drag coefficient [-]
        V_cruise               = 25,      # Cruise velocity [m/s]
    )
    
    print(Class_I_Weight_Estimate)
