import unittest
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from Weight_Estimations.Class_I_Weight_Estimation import Class_I_Weight_Estimation
from Weight_Estimations.Class_II_Weight_Estimation import ClassIIWeightEstimation


class TestClassIWeightFunctions(unittest.TestCase):

    def setUp(self):
        # Known inputs
        self.estimator = Class_I_Weight_Estimation(
            payload_weight_kg=500,
            residual_fuel_fraction=0.1,
            empty_weight_fraction=0.6,
            initial_mtow_guess_kg=5000,
            iteration_limit=100,
            tolerance=0.001,
            W1_WTO=0.99,
            W2_W1=0.98,
            W3_W2=0.99,
            W4_W3=0.995,
            W5_W4=0.985,
            W6_W5=0.99,
            W7_W6=0.99,
            W8_W7=0.995,
            Wfinal_W8=0.99,
            battery_power_available=10000,  # 10kW
            battery_specific_energy_Wh_per_kg=200,
            endurance=1,  # 1 day
            n_p=0.82,
            c_p=0.3,
            g=9.80665,
            A=20,
            e=0.85,
            CD0=0.025
        )

    def test_determine_mtow_convergence(self):
        self.assertTrue(self.estimator.converged, "MTOW should converge.")
        self.assertAlmostEqual(
            self.estimator.estimated_MTOM,
            self.estimator.results["Maximum Take-off Weight [kg]"],
            delta=0.01
        )

    def test_determine_used_fuel(self):
        M_ff = (
            0.99 * 0.98 * 0.99 * 0.995 * 0.985 * 0.99 * 0.99 * 0.995 * 0.99
        )
        expected_used_fuel = (1 - M_ff) * self.estimator.MTOW_kg
        self.estimator.Determine_Used_Fuel()
        self.assertAlmostEqual(
            self.estimator.W_f_used,
            expected_used_fuel,
            delta=0.01
        )

    def test_determine_mlw(self):
        self.estimator.Determine_MLW()
        oew = self.estimator.estimated_OEM
        payload = self.estimator.payload_weight_kg
        battery = self.estimator.battery_mass_kg
        residual_fuel = self.estimator.estimated_fuel_mass * 0.1
        expected_mlw = oew + payload + battery + residual_fuel
        self.assertAlmostEqual(
            self.estimator.results["Maximum Landing Weight [kg]"],
            round(expected_mlw, 3),
            delta=0.01
        )

    def test_battery_mass_calculation(self):
        expected_battery_mass = (10000 * 24) / 200
        self.assertAlmostEqual(
            self.estimator.battery_mass_kg,
            expected_battery_mass,
            delta=0.001
        )

    def test_determine_l_d_ratio(self):
        self.estimator.Determine_Maximum_Lift_Drag_Ratio()
        expected_l_d = ((3.1416 * 20 * 0.85) / (4 * 0.025)) ** 0.5
        self.assertAlmostEqual(
            self.estimator.L_D,
            expected_l_d,
            delta=0.01
        )

if __name__ == "__main__":
    unittest.main()
