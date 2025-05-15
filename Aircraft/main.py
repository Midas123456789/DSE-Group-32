# main.py

from Requirements import Requirements
from Atmosphere.ISA_Calculator import ISA_Calculator

from Weight_Estimations.Class_I_Weight_Estimation import Class_I_Weight_Estimation
from Weight_Estimations.Class_II_Weight_Estimation import ClassIIWeightEstimation

from Aerodynamics.Aircraft_Aerodynamics import AircraftAerodynamic
from Aerodynamics.Aero_plotting import plot_h_Preq, plot_A_LD, plot_feasible_S_V, plot_h_V, plot_h_W, plot_P_V

from Weight_Estimations.WP_WS_diagram import WP_WS_Diagram
from Performance.Performance import Performance


class AircraftInputs:
    """Container for all aircraft design input parameters."""
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)


class Aircraft:
    """Represents an aircraft with performance, aerodynamic and weight estimation subsystems."""
    def __init__(self, inputs: AircraftInputs):
        self.inputs = inputs

        # ISA Properties
        self.ISA = ISA_Calculator(
            altitude=inputs.h_cruise,
            length=inputs.chord_length,
        )
        self.g = self.ISA.results[inputs.h_cruise]["Gravity [m/s2]"]
        self.rho = self.ISA.results[inputs.h_cruise]["Density [kg/m3]"]

        # Class I weight estimation
        self.class_I = Class_I_Weight_Estimation(
            payload_weight_kg=inputs.payload_weight_kg,
            empty_weight_fraction=inputs.empty_weight_fraction,
            initial_mtow_guess_kg=inputs.initial_mtow_guess_kg,
            residual_fuel_fraction=inputs.residual_fuel_fraction,
            W1_WTO=inputs.W1_WTO, W2_W1=inputs.W2_W1, W3_W2=inputs.W3_W2,
            W4_W3=inputs.W4_W3, W5_W4=inputs.W5_W4, W6_W5=inputs.W6_W5,
            W7_W6=inputs.W7_W6, W8_W7=inputs.W8_W7, Wfinal_W8=inputs.Wfinal_W8,
            battery_power_available=inputs.battery_power_available,
            endurance=inputs.endurance,
            battery_specific_energy_Wh_per_kg=inputs.battery_specific_energy_Wh_per_kg,
            n_p=inputs.n_p, c_p=inputs.c_p,
            A=inputs.A, e=inputs.e, CD0=inputs.CD0, g=self.g,
            iteration_limit=100, tolerance=0.01,
        )

        self.aero = AircraftAerodynamic(
            W=self.class_I.estimated_MTOM,
            h=inputs.h_cruise, n_p=inputs.n_p,
            S=inputs.S, A=inputs.A, e=inputs.e,
            CD0=inputs.CD0, CL=inputs.CL_cruise
        )

        # Performance and aerodynamics
        self.performance = Performance(self)

if __name__ == "__main__":
    req = Requirements()
    inputs = AircraftInputs(
        
        # Mission
        endurance=req.endurance,
        
        # Weights
        payload_weight_kg=req.payload_weight_kg,
        initial_mtow_guess_kg=10 * req.payload_weight_kg,
        empty_weight_fraction=0.2,
        
        # Batteries
        battery_power_available=req.power_available,
        battery_specific_energy_Wh_per_kg=435,
        
        # Currently not using fuel!!!
        residual_fuel_fraction=0.00,
        W1_WTO=0.997,
        W2_W1=0.995,
        W3_W2=0.996,
        W4_W3=0.998,
        W5_W4=0.931,
        W6_W5=0.999,
        W7_W6=0.998,
        W8_W7=0.995,
        Wfinal_W8=0.997,
        
        # Propulsion system
        n_p=0.85, 
        c_p=0.6,
        sfc_kg_per_kw_hr=0.2,
        
        # Design parameters
        CL_cruise=1.2, 
        CL_TO=2.0, 
        CL_land=2.2,
        
        # Configuration
        CD0=0.020, 
        S=30, 
        A=25, 
        e=0.85,
        chord_length=3, 
        #fuselage_length=4, 
        #fuselage_diameter=0.001,
        
        # Choice
        h_cruise=15000,
        propulsion_type='battery',
        
    )
    ac = Aircraft(inputs)
    
    #print(50 * '-')
    #print(f"MTOM [kg]: {ac.class_I.estimated_MTOM:.2f} kg")
    #print(f"Payload mass [kg]: {req.payload_weight_kg:.2f} kg")
    #print("")
    #print(f"Battery mass required to fly for {req.endurance * 24} hours without recharging [kg]: {ac.performance.battery_mass_kg:.2f} kg")
    #print(f"Hydrogen fuel mass required to fly for {req.endurance * 24} hours without refuelling [kg]: {ac.performance.hydrogen_mass_kg:.2f} kg")
    #print(f"Hydrogen tank mass required to fit {ac.performance.hydrogen_mass_kg:.2f} kg: {ac.performance.hydrogen_tank_mass_kg:.2f} kg")
    #print(f"Cruise velocity for maximum endurance [m/s]: {ac.performance.optimum_V:.2f} m/s")
    #print(50 * '-')
    
    print(ac.ISA)
    print(ac.class_I)
    print(ac.aero)
    print(ac.performance)
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    """
    WP_WS_Diagram(V_cruise       = ac.performance.optimum_V, 
                V_stall          = 30, 
                climb_rate       = 5, 
                climb_gradient   = 0.01, 
                e                = ac.inputs.e, 
                CD0              = ac.inputs.CD0, 
                h                = ac.inputs.h_cruise, 
                CL_MAX_clean     = ac.inputs.CL_cruise, 
                CL_MAX_land      = ac.inputs.CL_land, 
                CL_MAX_TO        = ac.inputs.CL_TO, 
                A_design         = ac.inputs.A, 
                TOP              = 150, 
                landing_distance = 700, 
                f                = 1, 
                n_p              = ac.inputs.n_p, 
                prop_setting     = 1,
                rho              = ac.rho,
                )"""