import unittest
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from Weight_Estimations.Class_I_Weight_Estimation import Class_I_Weight_Estimation  # Replace with the actual module import
import numpy as np

class TestClassIWeightEstimation(unittest.TestCase):

    def setUp(self):
        # Set up basic parameters for testing
        self.payload_weight_kg = 1000
        self.residual_fuel_fraction = 0.0
        self.empty_weight_fraction = 0.7
        self.initial_mtow_guess_kg = 15000
        self.iteration_limit = 100
        self.tolerance = 1.0
        self.W1_WTO = 0.99
        self.W2_W1 = 0.99
        self.W3_W2 = 0.99
        self.W4_W3 = 0.99
        self.W5_W4 = 0.99
        self.W6_W5 = 0.99
        self.W7_W6 = 0.99
        self.W8_W7 = 0.99
        self.Wfinal_W8 = 0.99
        self.n_p = 0.8
        self.c_p = 0.3
        self.g = 9.80665
        self.battery_power_available = 5000
        self.battery_specific_energy_Wh_per_kg = 435
        self.endurance = 1
        self.charge_endurance = 1
        self.propulsion_type = 'hydrogen'

    def test_mtow_convergence(self):
        # Create an instance of Class_I_Weight_Estimation
        weight_estimator = Class_I_Weight_Estimation(
            payload_weight_kg=self.payload_weight_kg,
            residual_fuel_fraction=self.residual_fuel_fraction,
            empty_weight_fraction=self.empty_weight_fraction,
            initial_mtow_guess_kg=self.initial_mtow_guess_kg,
            iteration_limit=self.iteration_limit,
            tolerance=self.tolerance,
            W1_WTO=self.W1_WTO,
            W2_W1=self.W2_W1,
            W3_W2=self.W3_W2,
            W4_W3=self.W4_W3,
            W5_W4=self.W5_W4,
            W6_W5=self.W6_W5,
            W7_W6=self.W7_W6,
            W8_W7=self.W8_W7,
            Wfinal_W8=self.Wfinal_W8,
            n_p=self.n_p,
            c_p=self.c_p,
            g=self.g,
            battery_power_available=self.battery_power_available,
            battery_specific_energy_Wh_per_kg=self.battery_specific_energy_Wh_per_kg,
            endurance=self.endurance,
            charge_endurance=self.charge_endurance,
            propulsion_type=self.propulsion_type
        )

        # Run MTOW determination
        weight_estimator.Determine_MTOW()

        # Test if MTOW converged within the specified tolerance
        self.assertTrue(weight_estimator.converged)
        self.assertGreater(weight_estimator.estimated_MTOM, 0)

    def test_fuel_and_battery_calculation(self):
        # Create an instance of Class_I_Weight_Estimation
        weight_estimator = Class_I_Weight_Estimation(
            payload_weight_kg=self.payload_weight_kg,
            residual_fuel_fraction=self.residual_fuel_fraction,
            empty_weight_fraction=self.empty_weight_fraction,
            initial_mtow_guess_kg=self.initial_mtow_guess_kg,
            iteration_limit=self.iteration_limit,
            tolerance=self.tolerance,
            W1_WTO=self.W1_WTO,
            W2_W1=self.W2_W1,
            W3_W2=self.W3_W2,
            W4_W3=self.W4_W3,
            W5_W4=self.W5_W4,
            W6_W5=self.W6_W5,
            W7_W6=self.W7_W6,
            W8_W7=self.W8_W7,
            Wfinal_W8=self.Wfinal_W8,
            n_p=self.n_p,
            c_p=self.c_p,
            g=self.g,
            battery_power_available=self.battery_power_available,
            battery_specific_energy_Wh_per_kg=self.battery_specific_energy_Wh_per_kg,
            endurance=self.endurance,
            charge_endurance=self.charge_endurance,
            propulsion_type=self.propulsion_type
        )

        # Run MTOW determination and fuel estimation
        weight_estimator.Determine_MTOW()
        weight_estimator.Determine_Used_Fuel()

        # Test if fuel weight and battery mass are calculated
        self.assertGreater(weight_estimator.estimated_fuel_mass, 0)
        self.assertGreater(weight_estimator.battery_mass_kg, 0)

    def test_output_string_format(self):
        # Create an instance of Class_I_Weight_Estimation
        weight_estimator = Class_I_Weight_Estimation(
            payload_weight_kg=self.payload_weight_kg,
            residual_fuel_fraction=self.residual_fuel_fraction,
            empty_weight_fraction=self.empty_weight_fraction,
            initial_mtow_guess_kg=self.initial_mtow_guess_kg,
            iteration_limit=self.iteration_limit,
            tolerance=self.tolerance,
            W1_WTO=self.W1_WTO,
            W2_W1=self.W2_W1,
            W3_W2=self.W3_W2,
            W4_W3=self.W4_W3,
            W5_W4=self.W5_W4,
            W6_W5=self.W6_W5,
            W7_W6=self.W7_W6,
            W8_W7=self.W8_W7,
            Wfinal_W8=self.Wfinal_W8,
            n_p=self.n_p,
            c_p=self.c_p,
            g=self.g,
            battery_power_available=self.battery_power_available,
            battery_specific_energy_Wh_per_kg=self.battery_specific_energy_Wh_per_kg,
            endurance=self.endurance,
            charge_endurance=self.charge_endurance,
            propulsion_type=self.propulsion_type
        )

        # Run MTOW determination
        weight_estimator.Determine_MTOW()

        # Capture the string output
        result_str = str(weight_estimator)

        # Check if the string contains expected parameters
        self.assertIn("Class I Weight Estimation Results:", result_str)
        self.assertIn("Maximum Take-off Mass [kg]", result_str)
        self.assertIn("Operating Empty Mass [kg]", result_str)
        self.assertIn("Fuel Mass [kg]", result_str)
        self.assertIn("Battery Mass [kg]", result_str)


if __name__ == "__main__":
    unittest.main()
