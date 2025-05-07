import math

class ISA_Calculator:
    
    def __init__(self):
        
        # Constants
        self.g0    = 9.80665        # Sea level gravitational acceleration (m/s^2)
        self.R     = 287.05287      # Specific gas constant for dry air (J/kg·K)
        self.T0    = 288.15         # Sea level standard temperature (K)
        self.P0    = 101325.0       # Sea level standard pressure (Pa)
        self.rho0  = 1.225          # Sea level standard density (kg/m^3)
        
        # ISA temperature lapse rates (K/m) and base altitudes (m)
        layers = [
            (0,     11000, -0.0065),  # Troposphere
            (11000, 20000,  0.0),     # Tropopause
            (20000, 32000,  0.001),   # Lower Stratosphere
            (32000, 47000,  0.0028),  # Higher Stratosphere
            (47000, 51000,  0.0),     # Stratopause
            (51000, 71000, -0.0028),  # Lower Mesosphere
            (71000, 84852, -0.002),   # Higher Mesosphere
        ]
    
    def define_properties(h):
        if h < 0 or h > 84852:
            raise ValueError("Altitude must be between 0 and 84,852 meters.")
        
        T = T0
        P = P0
        for h_base, h_top, lapse in layers:
            if h < h_top:
                if lapse == 0:
                    T = T
                    P *= math.exp(-g0 / (R * T) * (h - h_base))
                else:
                    T = T + lapse * (h - h_base)
                    P *= (T / (T - lapse * (h - h_base))) ** (-g0 / (R * lapse))
                break
            else:
                if lapse == 0:
                    P *= math.exp(-g0 / (R * T) * (h_top - h_base))
                else:
                    T = T + lapse * (h_top - h_base)
                    P *= (T / (T - lapse * (h_top - h_base))) ** (-g0 / (R * lapse))
        
        rho = P / (R * T)
        g = g0 * (6371000 / (6371000 + h))**2
        return {
            "Temperature (K)": T,
            "Pressure (Pa)": P,
            "Density (kg/m³)": rho,
            "Gravity (m/s²)": g
        }

# Example usage
altitude = 20000
isa = ISA_Calculator()
for k, v in isa.define_properties(altitude).items():
    print(f"{k}: {v:.3f}")
