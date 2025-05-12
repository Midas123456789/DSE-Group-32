# main.py

from Class_I_Weight_Estimation import Class_I_Weight_Estimation
from Aircraft_Aerodynamics import AircraftAerodynamic
from Requirements import Requirements
from ISA_Calculator import ISA_Calculator
from Aero_plotting import plot_h_Preq, plot_A_LD, plot_feasible_S_V, plot_h_V, plot_h_W, plot_P_V
from Performance import Performance
from Class_II_Weight_Estimation import ClassIIWeightEstimation


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
            battery_specific_energy_Wh_per_kg=inputs.battery_specific_energy_Wh_per_kg,
            n_p=inputs.n_p, c_p=inputs.c_p,
            A=inputs.A, e=inputs.e, CD0=inputs.CD0, g=self.g
        )

        print(f"[0] MTOM: {self.class_I.estimated_MTOM:.2f} kg | "
              f"OEW_frac: {(self.class_I.estimated_OEM / self.class_I.estimated_MTOM):.3f} | "
              f"Battery mass: {self.class_I.battery_mass_kg:.2f} kg | Δ: None")

        # Performance and aerodynamics
        self.performance = Performance(self)
        self.aero = AircraftAerodynamic(
            W=self.class_I.estimated_MTOM,
            h=inputs.h_cruise,
            V=self.performance.optimum_V,
            S=inputs.S, A=inputs.A, e=inputs.e,
            CD0=inputs.CD0, CL=inputs.CL_cruise
        )

        # Class II weight estimation
        self.class_II = ClassIIWeightEstimation(
            payload_weight_kg=inputs.payload_weight_kg,
            battery_mass_kg=self.class_I.battery_mass_kg,
            S=inputs.S, A=inputs.A, g=self.g,
            fuselage_length=inputs.fuselage_length,
            fuselage_diameter=inputs.fuselage_diameter,
            initial_mtow_guess_kg=self.class_I.estimated_MTOM,
            W1_WTO=inputs.W1_WTO, W2_W1=inputs.W2_W1, W3_W2=inputs.W3_W2,
            W4_W3=inputs.W4_W3, W5_W4=inputs.W5_W4, W6_W5=inputs.W6_W5,
            W7_W6=inputs.W7_W6, W8_W7=inputs.W8_W7, Wfinal_W8=inputs.Wfinal_W8
        )

    def update_estimations(self, MTOW=None, OEW_fraction=None, battery_mass_kg=None):
        g = self.g

        initial_mtow_guess_kg = MTOW / g if MTOW else self.inputs.initial_mtow_guess_kg
        empty_weight_fraction = OEW_fraction if OEW_fraction else self.inputs.empty_weight_fraction

        if battery_mass_kg:
            battery_energy_Wh = battery_mass_kg * self.inputs.battery_specific_energy_Wh_per_kg
            battery_power_available = battery_energy_Wh
        else:
            battery_power_available = self.inputs.battery_power_available
            battery_mass_kg = battery_power_available / self.inputs.battery_specific_energy_Wh_per_kg

        self.class_I = Class_I_Weight_Estimation(
            payload_weight_kg=self.inputs.payload_weight_kg,
            empty_weight_fraction=empty_weight_fraction,
            initial_mtow_guess_kg=initial_mtow_guess_kg,
            residual_fuel_fraction=self.inputs.residual_fuel_fraction,
            W1_WTO=self.inputs.W1_WTO, W2_W1=self.inputs.W2_W1, W3_W2=self.inputs.W3_W2,
            W4_W3=self.inputs.W4_W3, W5_W4=self.inputs.W5_W4, W6_W5=self.inputs.W6_W5,
            W7_W6=self.inputs.W7_W6, W8_W7=self.inputs.W8_W7, Wfinal_W8=self.inputs.Wfinal_W8,
            battery_power_available=battery_power_available,
            battery_specific_energy_Wh_per_kg=self.inputs.battery_specific_energy_Wh_per_kg,
            n_p=self.inputs.n_p, c_p=self.inputs.c_p,
            A=self.inputs.A, e=self.inputs.e, CD0=self.inputs.CD0, g=g,
            battery_mass_kg=battery_mass_kg
        )

        self.performance = Performance(self)
        
        self.aero = AircraftAerodynamic(
            W=self.class_I.estimated_MTOM,
            h=self.inputs.h_cruise, V=self.performance.optimum_V,
            S=self.inputs.S, A=self.inputs.A, e=self.inputs.e,
            CD0=self.inputs.CD0, CL=self.inputs.CL_cruise
        )

        self.class_II = ClassIIWeightEstimation(
            payload_weight_kg=self.inputs.payload_weight_kg,
            battery_mass_kg=battery_mass_kg,
            S=self.inputs.S, A=self.inputs.A, g=g,
            fuselage_length=self.inputs.fuselage_length,
            fuselage_diameter=self.inputs.fuselage_diameter,
            initial_mtow_guess_kg=self.class_I.estimated_MTOM,
            W1_WTO=self.inputs.W1_WTO, W2_W1=self.inputs.W2_W1, W3_W2=self.inputs.W3_W2,
            W4_W3=self.inputs.W4_W3, W5_W4=self.inputs.W5_W4, W6_W5=self.inputs.W6_W5,
            W7_W6=self.inputs.W7_W6, W8_W7=self.inputs.W8_W7, Wfinal_W8=self.inputs.Wfinal_W8
        )


def weight_iteration(ac: Aircraft):
    """Iteratively refine weight and battery sizing until convergence."""
    convergence = float('inf')
    tolerance = 0.1
    tries = 0

    while convergence > tolerance:
        tries += 1

        previous_MTOM = ac.class_II.estimated_MTOM
        previous_OEW_fraction = ac.class_II.estimated_OEM / previous_MTOM

        drag = ac.aero.Drag()
        V = ac.performance.optimum_V
        power_required = drag * V
        energy_required_Wh = power_required * ac.performance.max_endurance / 3600
        battery_mass = energy_required_Wh / ac.inputs.battery_specific_energy_Wh_per_kg

        ac.update_estimations(MTOW=previous_MTOM, OEW_fraction=previous_OEW_fraction, battery_mass_kg=battery_mass)

        new_MTOM = ac.class_II.estimated_MTOM
        new_OEW_fraction = ac.class_II.estimated_OEM / new_MTOM
        convergence = abs(new_MTOM - previous_MTOM)

        print(f"[{tries}] MTOM: {new_MTOM:.2f} kg | OEW_frac: {new_OEW_fraction:.3f} | Battery mass: {battery_mass:.2f} kg | Δ: {convergence:.4f}")

        if tries > 100:
            print("Max iterations reached. Did not converge.")
            break


if __name__ == "__main__":
    req = Requirements()
    inputs = AircraftInputs(
        
        payload_weight_kg=req.payload_weight_kg,
        residual_fuel_fraction=0.00,
        empty_weight_fraction=0.7,
        initial_mtow_guess_kg=15 * req.payload_weight_kg,
        
        W1_WTO=1, 
        W2_W1=1, 
        W3_W2=1, 
        W4_W3=1, 
        W5_W4=1,
        W6_W5=1, 
        W7_W6=1, 
        W8_W7=1, 
        Wfinal_W8=1,
        
        n_p=0.85, 
        c_p=0.6,
        
        CL_cruise=1.2, 
        CL_TO=2.0, 
        CL_land=2.2,
        
        CD0=0.020, 
        
        chord_length=3, 
        S=30, 
        A=25, 
        e=0.85,
        
        h_cruise=15000,
        
        battery_power_available=20000,
        battery_specific_energy_Wh_per_kg=435,
        
        fuselage_length=2, 
        fuselage_diameter=0.1
        
    )

    ac = Aircraft(inputs)
    weight_iteration(ac)
    
    print(50 * '-')
    print(f"MTOM [kg]: {ac.class_II.estimated_MTOM:.2f} kg")
    print(f"Payload mass [kg]: {ac.class_II.payload_weight_kg:.2f} kg")
    print(f"Fuel mass used [kg]: {ac.class_I.W_f_used:.2f} kg")
    print(f"Battery mass [kg]: {ac.class_I.battery_mass_kg:.2f} kg")
    print(f"Cruise velocity [m/s]: {ac.performance.optimum_V:.2f} m/s")
    print(f"Endurance at cruise velocity [h]: {ac.performance.max_endurance / 3600:.2f} h")
    print(50 * '-')
