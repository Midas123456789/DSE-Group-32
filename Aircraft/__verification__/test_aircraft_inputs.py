import unittest
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from main import AircraftInputs

class TestAircraftInputs(unittest.TestCase):

    def test_aircraft_inputs_sets_attributes(self):
        # Create an instance of AircraftInputs
        inputs = AircraftInputs(chord    = 8, 
                                span     = 30,
                                sweep    = 28,
                                dihedral = 5,
                                )

        # Assert that the attributes are correctly set
        self.assertEqual(inputs.chord, 8)
        self.assertEqual(inputs.span, 30)
        self.assertEqual(inputs.sweep, 28)
        self.assertEqual(inputs.dihedral, 5)

if __name__ == "__main__":
    unittest.main()