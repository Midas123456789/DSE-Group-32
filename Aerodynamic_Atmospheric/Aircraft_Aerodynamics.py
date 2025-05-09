import matplotlib.pyplot as plt
import numpy as np
from Aerodynamics import *

class AircraftAerodynamic(Aerodynamic):
    def __init__(self, W=0, h=0, V=0, S=0, A=0, e=0, CD0=0, CL=0):
        super().__init__(weight=W, altitude=h)
        self.V = V
        self.S = S
        self.A = A 
        self.e = e 
        self.CD0 = CD0
        self.CL = CL

    @property
    def rho(self):
        return self.altitude_data[self.altitude]["Density [kg/m3]"]
    
    @property
    def a(self):
        return self.altitude_data[self.altitude]["Speed of Sound [m/s]"]
    
    @property
    def M(self):
        return self.V / self.a

    @property
    def CD(self):
        return self.drag_polar()

    def drag_polar(self, CL=None, CD0=None, A=None, e=None):
        CL = CL if CL is not None else self.CL
        CD0 = CD0 if CD0 is not None else self.CD0
        A = A if A is not None else self.A
        e = e if e is not None else self.e
        return (CD0 + (CL ** 2) / (np.pi * A * e))  # Drag Coefficient

    def Drag(self, CD=None, V=None, S=None, rho=None):
        CD = CD if CD is not None else self.CD
        V = V if V is not None else self.V
        S = S if S is not None else self.S
        rho = rho if rho is not None else self.rho
        return 0.5 * rho * (V ** 2) * S * CD

    def Lift(self, CL=None, V=None, S=None, rho=None):
        CL = CL if CL is not None else self.CL
        V = V if V is not None else self.V
        S = S if S is not None else self.S
        rho = rho if rho is not None else self.rho
        return 0.5 * rho * (V ** 2) * S * CL
    
    def min_velocity(self, CL=None, S=None, rho=None, weight=None):
        CL = CL if CL is not None else self.CL
        S = S if S is not None else self.S
        rho = rho if rho is not None else self.rho
        weight = weight if weight is not None else self.weight
        return np.sqrt((2 * weight) / (rho * S * CL))