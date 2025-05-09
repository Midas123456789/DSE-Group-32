from Class_I_Weight_Estimation import *
from Aircraft_Aerodynamics import *
from Requirements import *
from ISA_Calculator import *
from Aero_plotting import *

class Aircraft:
    
    def __init__(self, payload_weight_kg, residual_fuel_fraction ,empty_weight_fraction, initial_mtow_guess_kg, 
                W1_WTO, W2_W1, W3_W2, W4_W3, W5_W4, W6_W5, W7_W6, W8_W7, Wfinal_W8, n_p, c_p, A, e, C_D_0, V_cruise, 
                battery_energy_required_Wh, battery_specific_energy_Wh_per_kg, h_cruise, chord_length, S, CL_cruise, CL_land, CL_TO):
        
        # Requirements input
        self.payload_weight_kg = payload_weight_kg
        
        # Battery inputs
        self.battery_energy_required_Wh = battery_energy_required_Wh
        self.battery_specific_energy_Wh_per_kg = battery_specific_energy_Wh_per_kg
        
        # Inputs for first weight estimation
        self.residual_fuel_fraction = residual_fuel_fraction
        self.empty_weight_fraction = empty_weight_fraction
        self.initial_mtow_guess_kg = initial_mtow_guess_kg
        
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
        self.c_p = c_p
        self.C_D_0 = C_D_0
        self.V_cruise = V_cruise
        self.h_cruise = h_cruise
        
        # Coefficients
        self.CL_TO = CL_TO
        self.CL_cruise = CL_cruise
        self.CL_land = CL_land
        
        # Dimensions
        self.chord_length = chord_length
        self.S = S
        self.A = A
        self.e = e
        
        # Altitude plotting variables
        self.altitude_range = np.arange(0, 30000, 100).tolist()
        self.isa = ISA_Calculator(self.altitude_range)
        self.altitude_data = self.isa.results


if __name__ == "__main__":
    
    Requirements = Requirements()
    
    AC = Aircraft(
        
        payload_weight_kg  = Requirements.payload_weight_kg,
        
        residual_fuel_fraction  = 0.0, 
        empty_weight_fraction   = 0.7, 
        initial_mtow_guess_kg   = 15 * Requirements.payload_weight_kg, 
        
        W1_WTO     = 1, 
        W2_W1      = 1, 
        W3_W2      = 1, 
        W4_W3      = 1, 
        W5_W4      = 1, 
        W6_W5      = 1, 
        W7_W6      = 1, 
        W8_W7      = 1, 
        Wfinal_W8  = 1,
        
        n_p       = 0.85, 
        c_p       = 0.6, 
        
        CL_cruise = 1.5,
        CL_TO     = 2.0,
        CL_land   = 2.2,
        
        C_D_0     = 0.040, 
        V_cruise  = 25,
        h_cruise  = 20000,
        
        chord_length = 3,
        S            = 50,
        A         = 25, 
        e         = 0.8, 
        
        battery_energy_required_Wh         = 20000, 
        battery_specific_energy_Wh_per_kg  = 435,
        
    )
    
    ISA = ISA_Calculator(
        
        altitude  = AC.h_cruise,
        velocity  = AC.V_cruise,
        length    = AC.chord_length,
        
    )
    
    class_I = Class_I_Weight_Estimate = Class_I_Weight_Estimation(
        
        # Estimates and given weights
        payload_weight_kg      = AC.payload_weight_kg,
        empty_weight_fraction  = AC.empty_weight_fraction,
        initial_mtow_guess_kg  = AC.initial_mtow_guess_kg,
        
        # Fuel fractions
        residual_fuel_fraction = AC.residual_fuel_fraction,
        W1_WTO     = AC.W1_WTO,
        W2_W1      = AC.W2_W1,
        W3_W2      = AC.W3_W2,
        W4_W3      = AC.W4_W3,
        W5_W4      = AC.W5_W4,
        W6_W5      = AC.W6_W5,
        W7_W6      = AC.W7_W6,
        W8_W7      = AC.W8_W7,
        Wfinal_W8  = AC.Wfinal_W8,
        
        # Battery additions
        battery_energy_required_Wh        = AC.battery_energy_required_Wh, 
        battery_specific_energy_Wh_per_kg = AC.battery_specific_energy_Wh_per_kg,  
        
        # Properties for range and endurance
        n_p        = AC.n_p,
        c_p        = AC.c_p,
        A          = AC.A,
        e          = AC.e,
        C_D_0      = AC.C_D_0,
        V_cruise   = AC.V_cruise,
        g          = ISA.results[AC.h_cruise]["Gravity [m/s²]"],
    )
    
    Aerodynamics = AircraftAerodynamic(
        
        W    = class_I.estimated_MTOW, 
        h    = AC.h_cruise, 
        V    = AC.V_cruise, 
        S    = AC.S, 
        A    = AC.A, 
        e    = AC.e, 
        CD0  = AC.C_D_0, 
        CL   = AC.CL_cruise
    
    )
    
    