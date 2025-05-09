from Class_I_Weight_Estimation import Class_I_Weight_Estimation
from Aircraft_Aerodynamics import AircraftAerodynamic
from Requirements import Requirements
from ISA_Calculator import ISA_Calculator
from Aero_plotting import plot_h_Preq, plot_A_LD, plot_feasible_S_V, plot_h_V, plot_h_W, plot_P_V
from Performance import Performance

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
            length=inputs.chord_length,
        )
        self.rho = self.ISA.results[self.inputs.h_cruise]["Density [kg/m3]"]
        
        self.class_I = Class_I_Weight_Estimation(
            payload_weight_kg=inputs.payload_weight_kg,
            empty_weight_fraction=inputs.empty_weight_fraction,
            initial_mtow_guess_kg=inputs.initial_mtow_guess_kg,
            residual_fuel_fraction=inputs.residual_fuel_fraction,
            W1_WTO=inputs.W1_WTO, W2_W1=inputs.W2_W1, W3_W2=inputs.W3_W2,
            W4_W3=inputs.W4_W3, W5_W4=inputs.W5_W4, W6_W5=inputs.W6_W5,
            W7_W6=inputs.W7_W6, W8_W7=inputs.W8_W7, Wfinal_W8=inputs.Wfinal_W8,
            battery_power_available=inputs.battery_power_available,
            battery_specific_energy_Wh_per_kg=inputs.battery_specific_energy_Wh_per_kg,
            n_p=inputs.n_p, c_p=inputs.c_p,
            A=inputs.A, e=inputs.e,
            CD0=inputs.CD0,
            g=self.ISA.results[inputs.h_cruise]["Gravity [m/s2]"]
        )
        
        self.aero = AircraftAerodynamic(
            W=self.class_I.estimated_MTOW,
            h=inputs.h_cruise,
            S=inputs.S, A=inputs.A, e=inputs.e,
            CD0=inputs.CD0, CL=inputs.CL_cruise
        )
        
        self.performance = Performance(self)


if __name__ == "__main__":
    
    req = Requirements()
    
    inputs = Aircraft_Inputs(
        
        payload_weight_kg  = req.payload_weight_kg,
        
        residual_fuel_fraction  = 0.01, 
        empty_weight_fraction   = 0.7, 
        initial_mtow_guess_kg   = 15 * req.payload_weight_kg, 
        
        W1_WTO     = 0.997, 
        W2_W1      = 0.995, 
        W3_W2      = 0.996, 
        W4_W3      = 0.998, 
        W5_W4      = 0.931, 
        W6_W5      = 0.999, 
        W7_W6      = 0.998, 
        W8_W7      = 0.995, 
        Wfinal_W8  = 0.997,
        
        n_p       = 0.85, 
        c_p       = 0.3, 
        
        CL_cruise = 1.2,
        CL_TO     = 2.0,
        CL_land   = 2.2,
        
        CD0       = 0.020, 
        h_cruise  = 15000,
        
        chord_length = 3,
        S            = 30,
        A            = 25, 
        e            = 0.85, 
        
        battery_power_available            = 20000, 
        battery_specific_energy_Wh_per_kg  = 435,
    
    )
    
    ac = Aircraft(inputs)
    
    # Example usage TODO: Input found cruise speed into aerodynamics
    print(50*'-')
    print(f"MTOM according to Class I [kg]: {ac.class_I.estimated_MTOM:.2f} kg")
    print(f"Fuel mass used according to Class I [kg]: {ac.class_I.W_f_used:.2f} kg")
    print(f"Battery mass according to Class I [kg]: {ac.class_I.battery_mass_kg} kg")
    ac.performance.plot_endurance_vs_velocity()
    print(f"Optimum cruise velocity for maximum endurance [m/s]: {ac.performance.optimum_V:.2f} m/s")
    print(f"Maximum endurance at optimum velocity [s]: {ac.performance.max_endurance:.2f} s")
    print(50*'-')