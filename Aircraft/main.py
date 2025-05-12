# main.py

from Requirements import Requirements

from Weight_Estimations.Class_I_Weight_Estimation import Class_I_Weight_Estimation
from Weight_Estimations.Class_II_Weight_Estimation import ClassIIWeightEstimation

from Aerodynamics.Aircraft_Aerodynamics import AircraftAerodynamic
from Aerodynamics.Aero_plotting import plot_h_Preq, plot_A_LD, plot_feasible_S_V, plot_h_V, plot_h_W, plot_P_V

from ISA_Calculator import ISA_Calculator
from Others.Performance import Performance


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

        print(f"[0] MTOM: {self.class_I.estimated_MTOM:.2f} kg | "
            f"OEW_frac: {(self.class_I.estimated_OEM / self.class_I.estimated_MTOM):.3f} | "
            f"Battery mass: {self.class_I.battery_mass_kg:.2f} kg | Δ: None")

        self.aero = AircraftAerodynamic(
            W=self.class_I.estimated_MTOM,
            h=inputs.h_cruise,
            S=inputs.S, A=inputs.A, e=inputs.e,
            CD0=inputs.CD0, CL=inputs.CL_cruise
        )

        # Performance and aerodynamics
        self.performance = Performance(self)
        
        # Class II weight estimation
        self.class_II = ClassIIWeightEstimation(
            payload_weight_kg=inputs.payload_weight_kg,
            battery_mass_kg=self.performance.battery_mass_kg,
            S=inputs.S, A=inputs.A, g=self.g,
            fuselage_length=inputs.fuselage_length,
            fuselage_diameter=inputs.fuselage_diameter,
            initial_mtow_guess_kg=self.class_I.estimated_MTOM,
            W1_WTO=inputs.W1_WTO, W2_W1=inputs.W2_W1, W3_W2=inputs.W3_W2,
            W4_W3=inputs.W4_W3, W5_W4=inputs.W5_W4, W6_W5=inputs.W6_W5,
            W7_W6=inputs.W7_W6, W8_W7=inputs.W8_W7, Wfinal_W8=inputs.Wfinal_W8
        )

if __name__ == "__main__":
    req = Requirements()
    inputs = AircraftInputs(
        
        # Mission
        endurance=req.endurance,
        
        # Weights
        payload_weight_kg=req.payload_weight_kg,
        initial_mtow_guess_kg=15 * req.payload_weight_kg,
        empty_weight_fraction=0.7,
        
        # Batteries
        battery_power_available=req.power_available,
        battery_specific_energy_Wh_per_kg=435,
        
        # Currently not using fuel!!!
        residual_fuel_fraction=0.00,
        W1_WTO=1, W2_W1=1, W3_W2=1, W4_W3=1, W5_W4=1,W6_W5=1, W7_W6=1, W8_W7=1, Wfinal_W8=1,
        
        # Propulsion system
        n_p=0.85, 
        c_p=0.6,
        
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
        fuselage_length=2, 
        fuselage_diameter=0.001,
        
        # Choice
        h_cruise=15000,
        
    )
    ac = Aircraft(inputs)
    
    convergence_rate = 0.01
    break_value = 100
    counter = 0
    while counter <= break_value:
        counter += 1
    
    ac.performance.plot_drag_to_power_for_max_endurance()
    
    print(50 * '-')
    print(f"MTOM [kg]: {ac.class_II.estimated_MTOM:.2f} kg")
    print(f"Payload mass [kg]: {ac.class_II.payload_weight_kg:.2f} kg")
    print(f"Fuel mass used [kg]: {ac.class_II.estimated_fuel_weight_kg:.2f} kg")
    print(f"Battery mass [kg]: {ac.class_II.estimated_battery_mass_kg:.2f} kg")
    print("")
    print(f"Battery power required [W]: {ac.performance.battery_power_required:.2f} W")
    print(f"Cruise velocity for maximum endurance [m/s]: {ac.performance.optimum_V:.2f} m/s")
    print(50 * '-')
    
    # TODO: fix the weird values for P_req and create iteration
    # The weird values were caused by the fact you made it like you're flying 7 days on one charge!!! you should differentiate this
    # Implement recharging
    