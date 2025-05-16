import aerosandbox as asb

from Requirements import Requirements
from Weight_Estimations.Class_I_Weight_Estimation import Class_I_Weight_Estimation
from Performance.Performance import Performance


class AircraftInputs:
    """Container for all aircraft design input parameters."""
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)


class Aircraft:
    """Represents an aircraft with performance, aerodynamic and weight estimation subsystems."""
    def __init__(self, inputs: AircraftInputs):
        
        # Set-up toolbox tools
        self.opti = asb.Opti()
        self.sol = self.opti.solve(verbose=False)
        
        # Get input parameters
        self.inputs = inputs
        
        # Run atmospheric model
        self.ISA = asb.Atmosphere(altitude=inputs.h_cruise)
        self.rho = self.ISA.density()
        print(f"Atmospheric model tuned at {inputs.h_cruise}")
        
        # Set up aerodynamic profile
        wing_airfoil = asb.Airfoil("sd7037")
        wing = asb.Wing(
            name       = "Main",
            symmetric  = True,
            xsecs=[
                asb.WingXSec(  # 0
                    xyz_le=[0, 0, 0],
                    chord=inputs.chord,
                    airfoil=wing_airfoil,
                ),
                asb.WingXSec(  # 1
                    xyz_le=[inputs.sweep, inputs.span * (3/4), inputs.dihedral],
                    chord=inputs.chord,
                    airfoil=wing_airfoil,
                ),
                asb.WingXSec(  # 2
                    xyz_le=[inputs.sweep * (7/8), inputs.span, inputs.dihedral],
                    chord=inputs.chord,
                    airfoil=wing_airfoil,
                )
            ]
        )
        self.configuration = asb.Airplane(wings=[wing], fuselages=[])
        print(f"Wing & Fuselage configuration is defined")

        # Class I weight estimation
        print(f"Performing Weight Estimation model")
        self.class_I = Class_I_Weight_Estimation(
            payload_weight_kg      = inputs.payload_weight_kg,
            empty_weight_fraction  = inputs.empty_weight_fraction,
            initial_mtow_guess_kg  = inputs.initial_mtow_guess_kg,
            residual_fuel_fraction = inputs.residual_fuel_fraction,
            
            W1_WTO = inputs.W1_WTO,  W2_W1 = inputs.W2_W1, W3_W2     = inputs.W3_W2,
            W4_W3  = inputs.W4_W3,   W5_W4 = inputs.W5_W4, W6_W5     = inputs.W6_W5,
            W7_W6  = inputs.W7_W6,   W8_W7 = inputs.W8_W7, Wfinal_W8 = inputs.Wfinal_W8,
            
            battery_specific_energy_Wh_per_kg  = inputs.battery_specific_energy_Wh_per_kg,
            battery_power_available  = inputs.battery_power_available,
            endurance  = inputs.endurance,
            
            n_p = inputs.n_p, 
            c_p = inputs.c_p,
            g   = inputs.g,
            
            iteration_limit = 100, 
            tolerance       = 0.01,
        )
        
        # Performance
        print(f"Starting Performance model")
        self.performance = Performance(ac=self, plot=True)
        print(f"Ending Performance model")


if __name__ == "__main__":
    req = Requirements()
    
    inputs = AircraftInputs(
        
        # Mission
        endurance=req.endurance,
        g = 9.80665,
        
        # Weights
        payload_weight_kg=req.payload_weight_kg,
        initial_mtow_guess_kg=10 * req.payload_weight_kg,
        empty_weight_fraction=0.2,
        
        # Batteries
        battery_power_available=req.power_available,
        battery_specific_energy_Wh_per_kg=435,
        
        # Hydrogen
        hydrogen_specific_energy_Wh_per_kg = 33333, # 120 MJ/kg in Wh
        hydrogen_density_kg_per_m3 = 42,
        tank_mass_fraction = 0.15,
        
        # Currently not using fuel!!!
        residual_fuel_fraction=0.00,
        W1_WTO    = 0.997,
        W2_W1     = 0.995,
        W3_W2     = 0.996,
        W4_W3     = 0.998,
        W5_W4     = 0.931,
        W6_W5     = 0.999,
        W7_W6     = 0.998,
        W8_W7     = 0.995,
        Wfinal_W8 = 0.997,
        
        # Propulsion system
        n_p=0.85, 
        c_p=0.6,
        
        # Design parameters
        CL_cruise = 1.2, 
        CL_TO     = 2.0, 
        CL_land   = 2.2,
        
        # Configuration
        chord    = 8, 
        span     = 30,
        sweep    = 28,
        dihedral = 5,
        
        # Choice
        h_cruise        = 15000,
        propulsion_type = 'battery',
        
    )
    ac = Aircraft(inputs)
    
    print(ac.class_I)
    print(ac.performance)
    