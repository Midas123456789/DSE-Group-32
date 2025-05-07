import math

class ISA_Calculator:
    
    def __init__(self, altitude):
        
        # Constants
        self.g0    = 9.80665         # Sea level gravitational acceleration [m/s^2]
        self.R     = 287.05287       # Specific gas constant for dry air [J/kg·K]
        self.T0    = 288.15          # Sea level standard temperature [K]
        self.P0    = 101325.0        # Sea level standard pressure [Pa]
        self.rho0  = 1.225           # Sea level standard density [kg/m^3]
        self.earth_radius = 6371000  # Radius of Earth [m]
        
        # ISA temperature alpha rates (K/m) and base altitudes (m)
        self.layers = [
            (0,     11000, -0.0065),  # Troposphere
            (11000, 20000,  0.0),     # Tropopause
            (20000, 32000,  0.001),   # Lower Stratosphere
            (32000, 47000,  0.0028),  # Higher Stratosphere
            (47000, 51000,  0.0),     # Stratopause
            (51000, 71000, -0.0028),  # Lower Mesosphere
            (71000, 84852, -0.002),   # Higher Mesosphere
        ]
        
        self.results = {}
        self.define_properties(altitude)
    
    def define_properties(self, altitude):
        if altitude < 0 or altitude > 84852:
            raise ValueError("Altitude must be between 0 and 84,852 meters.")
        
        T = self.T0
        P = self.P0
        for h_base, h_next, alpha in self.layers:
            if altitude < h_next:
                if alpha == 0: # Pause layers
                    T = T
                    P *= math.exp((-self.g0 / (self.R * T)) * (altitude - h_base))
                else: # Sphere layers
                    T = T + (alpha * (altitude - h_base))
                    P *= (T / (T - (alpha * (altitude - h_base)))) ** (-self.g0 / (self.R * alpha))
                break
            else:
                if alpha == 0: # Pause layers
                    P *= math.exp(-self.g0 / (self.R * T) * (h_next - h_base))
                else: # Sphere layers
                    T = T + alpha * (h_next - h_base)
                    P *= (T / (T - alpha * (h_next - h_base))) ** (-self.g0 / (self.R * alpha))
        
        # Define density and gravity at desired altitude
        rho = P / (self.R * T)
        g = self.g0 * (self.earth_radius / (self.earth_radius + altitude))**2
        self.results = {
            "Altitude [m]": altitude,
            "Temperature [K]": T,
            "Pressure [Pa]": P,
            "Density [kg/m³]": rho,
            "Gravity [m/s²]": g
        }
    
    def __str__(self):
        output = ["ISA Calculator Results:\n"]
        max_key_length = max(len(key) for key in self.results.keys())
        header = f"{'Parameter'.ljust(max_key_length)} | Value"
        output.append(header)
        output.append("-" * len(header))
        
        for key, value in self.results.items():
            output.append(f"{key.ljust(max_key_length)} | {value:>13,.2f}")
        return "\n".join(output)

# Example usage
isa = ISA_Calculator(altitude=20000)
print(isa)

