import math
import numpy as np

class ISA_Calculator:
    
    def __init__(self, altitude=0, velocity=0, length=0):
        # Constants
        self.g0 = 9.80665
        self.R = 287.05287
        self.T0 = 288.15
        self.P0 = 101325.0
        self.rho0 = 1.225
        self.earth_radius = 6371000
        self.gamma = 1.4
        self.C = 110.4
        self.mu_0 = 1.716e-5
        
        self.layers = [
            (0,     11000, -0.0065),
            (11000, 20000,  0.0),
            (20000, 32000,  0.001),
            (32000, 47000,  0.0028),
            (47000, 51000,  0.0),
            (51000, 71000, -0.0028),
            (71000, 84852, -0.002),
        ]
        
        self.results = {}
        
        if not isinstance(altitude, list):
            altitude = [altitude]
        for h in altitude:
            self.results[h] = {}  # Use altitude as key
            self.ISA_properties(h, self.results[h])
            self.mach_number(velocity, self.results[h])
            self.reynolds_number(velocity, length, self.results[h])
    
    def ISA_properties(self, altitude, result):
        if altitude < 0 or altitude > 84852:
            raise ValueError("Altitude must be between 0 and 84,852 meters.")
        
        T = self.T0
        P = self.P0
        
        for h_base, h_next, alpha in self.layers:
            if altitude < h_next:
                if alpha == 0:
                    P = self.pause_layer_calculation(T, P, altitude, h_base)
                else:
                    T, P = self.sphere_layer_calculation(T, P, altitude, h_base, alpha)
                break
            else:
                if alpha == 0:
                    P = self.pause_layer_calculation(T, P, h_next, h_base)
                else:
                    T, P = self.sphere_layer_calculation(T, P, h_next, h_base, alpha)
        
        rho = P / (self.R * T)
        g = self.g0 * (self.earth_radius / (self.earth_radius + altitude))**2
        
        result.update({
            "Temperature [K]": T,
            "Pressure [Pa]": P,
            "Density [kg/m3]": rho,
            "Gravity [m/s2]": g,
        })
    
    def sphere_layer_calculation(self, T, P, altitude, h_0, alpha):
        T = T + alpha * (altitude - h_0)
        P *= (T / (T - (alpha * (altitude - h_0)))) ** (-self.g0 / (self.R * alpha))
        return T, P
    
    def pause_layer_calculation(self, T, P, altitude, h_0):
        P *= math.exp((-self.g0 / (self.R * T)) * (altitude - h_0))
        return P
    
    def dynamic_viscosity(self, T):
        return self.mu_0 * (((T / self.T0)**1.5) * ((self.T0 + self.C) / (T + self.C)))
    
    def mach_number(self, V, result):
        a = math.sqrt(self.gamma * self.R * result["Temperature [K]"])
        result["Speed of Sound [m/s]"] = a
        result["Mach Number"] = V / a
    
    def reynolds_number(self, V, L, result):
        mu = self.dynamic_viscosity(result["Temperature [K]"])
        result["Reynolds Number"] = (result["Density [kg/m3]"] * V * L) / mu
    
    def __str__(self):
        output = ["ISA Calculator Results:\n"]
        for altitude, res in self.results.items():
            output.append(f"Altitude: {altitude} m")
            max_key_length = max(len(key) for key in res)
            header = f"{'Parameter'.ljust(max_key_length)} | Value"
            output.append(header)
            output.append("-" * len(header))
            for key, value in res.items():
                output.append(f"{key.ljust(max_key_length)} | {value:>13,.5f}")
            output.append("")
        return "\n".join(output)
    
    def plot(self):
        """ Visualize the ISA properties """
        pass


# Example Usage
if __name__ == "__main__":
    isa = ISA_Calculator(altitude=[15000, 30000, 50000], velocity=25, length=2)
    
    # Access temperature at 30,000 m:
    print(isa.results[15000]["Gravity [m/s2]"])
