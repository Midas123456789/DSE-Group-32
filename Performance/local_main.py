import aerosandbox as asb
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from Requirements import Requirements
from Weight_Estimations.Class_I_Weight_Estimation import Class_I_Weight_Estimation
from Performance import Performance
from Aerosandbox_models.LH2_tandem import *


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
        
        # Class I weight estimation
        self.class_I = Class_I_Weight_Estimation(
            
            propulsion_type = inputs.propulsion_type,
            
            payload_weight_kg      = inputs.payload_weight_kg,
            empty_weight_fraction  = inputs.empty_weight_fraction,
            initial_mtow_guess_kg  = inputs.initial_mtow_guess_kg,
            residual_fuel_fraction = inputs.residual_fuel_fraction,
            
            W1_WTO = inputs.W1_WTO,   W2_W1 = inputs.W2_W1,   W3_W2     = inputs.W3_W2,
            W4_W3  = inputs.W4_W3,    W5_W4 = inputs.W5_W4,   W6_W5     = inputs.W6_W5,
            W7_W6  = inputs.W7_W6,    W8_W7 = inputs.W8_W7,   Wfinal_W8 = inputs.Wfinal_W8,
            
            battery_specific_energy_Wh_per_kg  = inputs.battery_specific_energy_Wh_per_kg,
            battery_power_available  = inputs.battery_power_available,
            endurance        = inputs.endurance,
            charge_endurance = inputs.charge_endurance,
            
            n_p = inputs.n_p, 
            c_p = inputs.c_p,
            g   = inputs.g,
            
            iteration_limit = 100, 
            tolerance       = 0.01,
        )
        
        # Performance
        self.performance = Performance(
            
            plot = True,
            n_p = inputs.n_p,
            h_cruise = inputs.h_cruise,
            rho = self.rho,
            g = inputs.g,
            battery_power_available = inputs.battery_power_available,
            battery_specific_energy_Wh_per_kg = inputs.battery_specific_energy_Wh_per_kg,
            propulsion_type = inputs.propulsion_type,
            endurance = inputs.endurance,
            charge_endurance = inputs.charge_endurance,
            configuration = inputs.configuration,
            hydrogen_specific_energy_Wh_per_kg = inputs.hydrogen_specific_energy_Wh_per_kg, 
            hydrogen_density_kg_per_m3 = inputs.hydrogen_density_kg_per_m3, 
            tank_mass_fraction = inputs.tank_mass_fraction,
            S = inputs.S,
        )


if __name__ == "__main__":
    req = Requirements()
    
    span = 32.12
    N_cords = 2
    cords = [2.68, 0.7]
    wing_airfoil = asb.Airfoil("sd7037")
    
    y_sections = np.linspace(0, span / 2, N_cords)
    wing1 = asb.Wing(
        symmetric=True,
        xsecs=[
            asb.WingXSec(
                xyz_le=[-0.25 * cords[i], y_sections[i], 0],
                chord=cords[i],
                airfoil=wing_airfoil
            )
            for i in range(N_cords)
        ]
    )
    
    wing2 = asb.Wing(
        symmetric=True,
        xsecs=[
            asb.WingXSec(
                xyz_le=[-0.25 * cords[i], y_sections[i], 0],
                chord=cords[i],
                airfoil=wing_airfoil
            )
            for i in range(N_cords)
        ]
    ).translate([9,0,0])
    
    V_stab = [asb.Wing(
        name="V-stab",
        xsecs=[
            asb.WingXSec(
                xyz_le=[0, 0, 0],
                chord=3,
                airfoil=asb.Airfoil("ht08")
            ),
            asb.WingXSec(
                xyz_le=[2, 0, 2],
                chord=1,
                airfoil=asb.Airfoil("ht08")
            )
        ]
    ).translate([6, 4, 0]),

    asb.Wing(
        name="V-stab",
        xsecs=[
            asb.WingXSec(
                xyz_le=[0, 0, 0],
                chord=3,
                airfoil=asb.Airfoil("ht08")
            ),
            asb.WingXSec(
                xyz_le=[2, 0, 2],
                chord=1,
                airfoil=asb.Airfoil("ht08")
            )
        ]
    ).translate([6, -4, 0])
    ]
    
    fuselages=[
        asb.Fuselage(
            name="Fuselage",
            xsecs=[
                asb.FuselageXSec(
                    xyz_c=[xi * 11, 0, 0],
                    radius=5*asb.Airfoil("naca0020").local_thickness(x_over_c=xi)
                )
                for xi in np.cosspace(0, 1, 30)
            ]
        ).translate([-1,4,0]),
        
        asb.Fuselage(
            name="Fuselage",
            xsecs=[
                asb.FuselageXSec(
                    xyz_c=[xi * 11, 0, 0],
                    radius=5*asb.Airfoil("naca0020").local_thickness(x_over_c=xi)
                )
                for xi in np.cosspace(0, 1, 30)
            ]
        ).translate([-1,-4,0]),
        
        asb.Fuselage(
            name="Fuselage",
            xsecs=[
                asb.FuselageXSec(
                    xyz_c=[xi * 5, 0, 0],
                    radius=3*asb.Airfoil("naca0010").local_thickness(x_over_c=xi)
                )
                for xi in np.cosspace(0, 1, 30)
            ]
        ).translate([-1.5,0,-0.5])
    ]
    
    configuration = asb.Airplane(wings=[wing1] + [wing2] + V_stab, fuselages=fuselages)
    S = wing1.area() + wing2.area()

    inputs = AircraftInputs(
        
        # Mission
        endurance=req.endurance,
        charge_endurance=0.5,
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
        
        # Choice
        h_cruise        = 17000,
        propulsion_type = 'hydrogen', # battery (not fully integrated yet!)
        
        configuration = configuration,
        S = S,
    )
    ac = Aircraft(inputs)
    
    # Run functions
    ac.class_I.Determine_MTOW()
    ac.class_I.Determine_Used_Fuel()
    print(ac.class_I)
    
    ac.performance.optimize_for_maximum_endurance(ac.class_I.estimated_MTOM)
    #ac.performance.find_battery_or_hydrogen_weight_parameters()
    print(ac.performance)

    configuration.draw(backend="matplotlib",set_axis_visibility=False, ax=None, thin_wings=True)


