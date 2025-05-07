import numpy as np
import matplotlib.pyplot as plt
from ISA_Calculator import ISA_Calculator

class Aerodynamic:
    def __init__(self, weight, altitude):
        """
        Base aerodynamic model.
        Parameters:
        - weight: in Newtons
        - altitude: in meters
        """
        self.weight = weight
        self.altitude = altitude

        self.altitude_range = np.linspace(0, 30000, 100).tolist()
        self.isa = ISA_Calculator(self.altitude_range)
        self.altitude_data = self.isa.results
