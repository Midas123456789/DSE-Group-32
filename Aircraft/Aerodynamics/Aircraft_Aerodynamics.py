import matplotlib.pyplot as plt
import numpy as np
from Aerodynamics.Aerodynamics import Aerodynamic

class AircraftAerodynamic(Aerodynamic):
    def __init__(self, W=0, h=0, V=0, S=0, A=0, e=0, CD0=0, CL=0, n_p=1):
        super().__init__(weight=W, altitude=h)
        self.V = V
        self.S = S
        self.A = A 
        self.e = e 
        self.CD0 = CD0
        self.CL = CL
        self.n_p = n_p

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
    
    def compute_drag_and_power(self, V):
        """Returns total drag and power required at velocity V."""
        q = 0.5 * self.rho * V**2
        parasite_drag = q * self.S * self.CD0
        induced_drag = (self.weight**2) / (q * self.S * np.pi * self.A * self.e) if q > 0 else 0
        total_drag = parasite_drag + induced_drag
        power_required = total_drag * V / self.n_p  # in Watts
        return total_drag, power_required

    def __str__(self):
        output = ["\n"]
        output.append("Aircraft Aerodynamic Properties:")
        data = {
            "Altitude [m]": self.altitude,
            "Speed [m/s]": self.V,
            "Wing Area [m²]": self.S,
            "Aspect Ratio [-]": self.A,
            "Oswald Efficiency [-]": self.e,
            "CD0 [-]": self.CD0,
            "CL [-]": self.CL,
            "Mach Number [-]": self.M,
            "Air Density [kg/m³]": self.rho,
            "Speed of Sound [m/s]": self.a,
            "Drag Coefficient (CD)": self.CD,
            "Lift [N]": self.Lift(),
            "Drag [N]": self.Drag(),
            "Power Required [W]": self.compute_drag_and_power(self.V)[1],
        }

        max_key_len = max(len(k) for k in data)
        header = f"{'Parameter'.ljust(max_key_len)} | Value"
        output.append(header)
        output.append("-" * len(header))

        for key, val in data.items():
            try:
                output.append(f"{key.ljust(max_key_len)} | {val:>13,.3f}")
            except Exception:
                output.append(f"{key.ljust(max_key_len)} | {'N/A':>13}")

        return "\n".join(output)
