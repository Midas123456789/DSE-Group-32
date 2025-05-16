import unittest
import sys
import os
import math

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from Atmosphere.ISA_Calculator import ISA_Calculator


class TestISACalculator(unittest.TestCase):
    
    def setUp(self):
        self.calc = ISA_Calculator()  # No need to trigger calculations
        self.default_T = 288.15
        self.default_P = 101325.0
    
    def test_sphere_layer_calculation(self):
        T, P = self.calc.sphere_layer_calculation(T=self.default_T, P=self.default_P, altitude=11000, h_0=0, alpha=-0.0065)
        self.assertAlmostEqual(T, 216.65, delta=0.01)
        self.assertAlmostEqual(P, 22632.04, delta=0.01)
    
    def test_pause_layer_calculation(self):
        T = 216.65
        P = 22632.04
        P_out = self.calc.pause_layer_calculation(T=T, P=P, altitude=15000, h_0=11000)
        self.assertAlmostEqual(P_out, 12044.55, delta=0.01)
    
    def test_dynamic_viscosity(self):
        T = 250  # Typical temp at 10 km
        mu = self.calc.dynamic_viscosity(T)
        self.assertAlmostEqual(mu, 1.46e-5, delta=1e-6)
    
    def test_mach_number(self):
        result = {"Temperature [K]": 288.15}
        self.calc.mach_number(V=100, result=result)
        self.assertIn("Mach Number", result)
        self.assertAlmostEqual(result["Mach Number"], 100 / math.sqrt(self.calc.gamma * self.calc.R * 288.15), places=5)
    
    def test_reynolds_number(self):
        result = {
            "Temperature [K]": 288.15,
            "Density [kg/m3]": 1.225
        }
        self.calc.reynolds_number(V=50, L=2, result=result)
        self.assertIn("Reynolds Number", result)
        self.assertAlmostEqual(result["Reynolds Number"], 7138694.74431, delta=1)
    
    def test_gravity_variation(self):
        result = {}
        self.calc.ISA_properties(altitude=50000, result=result)
        g = result["Gravity [m/s2]"]
        self.assertAlmostEqual(g, 9.65452, delta=1e-5)


if __name__ == "__main__":
    unittest.main()
