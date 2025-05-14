import unittest
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from Weight_Estimations.Class_I_Weight_Estimation import Class_I_Weight_Estimation
from Weight_Estimations.Class_II_Weight_Estimation import ClassIIWeightEstimation
from Requirements import Requirements


class TestClassIWeightEstimation(unittest.TestCase):

    def setUp(self):
        self.req = Requirements()
        self.estimator = Class_I_Weight_Estimation(
            payload_weight_kg=self.req.payload_weight_kg,
            residual_fuel_fraction=0.0,
            empty_weight_fraction=0.7,
            initial_mtow_guess_kg=15 * self.req.payload_weight_kg,
            iteration_limit=100,
            tolerance=0.01,
            W1_WTO=1,
            W2_W1=1,
            W3_W2=1,
            W4_W3=1,
            W5_W4=1,
            W6_W5=1,
            W7_W6=1,
            W8_W7=1,
            Wfinal_W8=1,
            battery_power_available=20000,
            battery_specific_energy_Wh_per_kg=435,
            endurance=2,  # days
            n_p=0.82,
            c_p=0.3,
            g=9.80665,
            A=25,
            e=0.9,
            CD0=0.020
        )

    def test_convergence(self):
        self.assertTrue(self.estimator.converged, "MTOW iteration did not converge.")

    def test_battery_mass(self):
        expected = (20000 * 48) / 435
        self.assertAlmostEqual(self.estimator.battery_mass_kg, expected, delta=0.1)

    def test_used_fuel_estimate(self):
        self.assertAlmostEqual(self.estimator.W_f_used, 0.0, delta=1e-3)

    def test_mlw_greater_than_payload(self):
        mlw = self.estimator.results.get("Maximum Landing Weight [kg]", 0)
        self.assertGreater(mlw, self.req.payload_weight_kg)

    def test_mtow_greater_than_payload_plus_empty(self):
        mtow = self.estimator.results.get("Maximum Take-off Weight [kg]", 0)
        oew = self.estimator.results.get("Operating Empty Weight [kg]", 0)
        self.assertGreater(mtow, self.req.payload_weight_kg + oew)


class DummyClassI:
    """Mock class to simulate output from Class_I_Weight_Estimation"""
    battery_mass_kg = 120.0
    payload_weight_kg = 500.0
    MTOW_kg = 5000.0


class TestClassIIWeightEstimation(unittest.TestCase):

    def setUp(self):
        self.estimator = ClassIIWeightEstimation(
            payload_weight_kg=DummyClassI.payload_weight_kg,
            battery_mass_kg=DummyClassI.battery_mass_kg,
            initial_mtow_guess_kg=DummyClassI.MTOW_kg,
            fuselage_length=10,
            fuselage_diameter=1.5,
            S=30,
            A=25
        )

    def test_convergence(self):
        self.assertTrue(self.estimator.converged, "MTOW iteration did not converge.")

    def test_mtow_positive(self):
        self.assertGreater(self.estimator.MTOW_kg, 0, "MTOW should be greater than 0.")

    def test_component_keys_exist(self):
        keys = self.estimator.results.keys()
        for key in ["Wing [kg]", "Fuselage [kg]", "Empennage [kg]", "Landing Gear [kg]", "Avionics [kg]", "Propulsion [kg]"]:
            self.assertIn(key, keys)

    def test_oew_fraction_valid(self):
        self.assertGreater(self.estimator.estimated_OEM_fraction, 0.1)
        self.assertLess(self.estimator.estimated_OEM_fraction, 1.0)

    def test_total_mass_balance(self):
        mtow = self.estimator.estimated_MTOM
        sum_parts = (self.estimator.results["OEW [kg]"]
                    + self.estimator.results["Payload [kg]"]
                    + self.estimator.results["Battery [kg]"]
                    + self.estimator.results["Fuel [kg]"])
        self.assertAlmostEqual(mtow, sum_parts, delta=1.0)

if __name__ == "__main__":
    unittest.main()
