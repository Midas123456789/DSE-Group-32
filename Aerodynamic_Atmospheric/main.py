from Class_I_Weight_Estimation import *
from Aircraft_Aerodynamics import *
from Requirements import *
from ISA_Calculator import *
from Aero_plotting import *

class Aircraft_Inputs:

    def __init__(self, **kwargs):
        # Save all inputs as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)

class Aircraft:
    def __init__(self, inputs: Aircraft_Inputs):
        self.inputs = inputs
        
        # Subsystems initialized here
        self.ISA = ISA_Calculator(
            altitude=inputs.h_cruise,
            velocity=inputs.V_cruise,
            length=inputs.chord_length,
        )
        
        self.class_I = Class_I_Weight_Estimation(
            payload_weight_kg=inputs.payload_weight_kg,
            empty_weight_fraction=inputs.empty_weight_fraction,
            initial_mtow_guess_kg=inputs.initial_mtow_guess_kg,
            residual_fuel_fraction=inputs.residual_fuel_fraction,
            W1_WTO=inputs.W1_WTO, W2_W1=inputs.W2_W1, W3_W2=inputs.W3_W2,
            W4_W3=inputs.W4_W3, W5_W4=inputs.W5_W4, W6_W5=inputs.W6_W5,
            W7_W6=inputs.W7_W6, W8_W7=inputs.W8_W7, Wfinal_W8=inputs.Wfinal_W8,
            battery_energy_required_Wh=inputs.battery_energy_required_Wh,
            battery_specific_energy_Wh_per_kg=inputs.battery_specific_energy_Wh_per_kg,
            n_p=inputs.n_p, c_p=inputs.c_p,
            A=inputs.A, e=inputs.e,
            C_D_0=inputs.C_D_0, V_cruise=inputs.V_cruise,
            g=self.ISA.results[inputs.h_cruise]["Gravity [m/s²]"]
        )
        
        self.aero = AircraftAerodynamic(
            W=self.class_I.estimated_MTOW,
            h=inputs.h_cruise, V=inputs.V_cruise,
            S=inputs.S, A=inputs.A, e=inputs.e,
            CD0=inputs.C_D_0, CL=inputs.CL_cruise
        )


if __name__ == "__main__":
    
    req = Requirements()
    
    inputs = Aircraft_Inputs(
        
        payload_weight_kg  = req.payload_weight_kg,
        
        residual_fuel_fraction  = 0.0, 
        empty_weight_fraction   = 0.7, 
        initial_mtow_guess_kg   = 15 * req.payload_weight_kg, 
        
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
        
        CL_cruise = 1.2,
        CL_TO     = 2.0,
        CL_land   = 2.2,
        
        C_D_0     = 0.040, 
        V_cruise  = 25,
        h_cruise  = 20000,
        
        chord_length = 3,
        S            = 30,
        A         = 20, 
        e         = 0.7, 
        
        battery_energy_required_Wh         = 20000, 
        battery_specific_energy_Wh_per_kg  = 435,
        
    )
    
    ac = Aircraft(inputs)
    
    # Example usage
    print("MTOW:", ac.class_I.estimated_MTOW)
    print("Cruise Drag:", ac.aero.Drag())
    
    
    