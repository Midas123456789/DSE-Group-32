import math
import numpy as np

class ISA_Calculator:
    """
    Computes ISA atmospheric properties and derived aerodynamic parameters (Mach, Reynolds).
    
    Args:
        altitude (float or list of float): Altitude(s) in meters.
        velocity (float): Velocity in m/s for Mach and Reynolds calculations.
        length (float): Reference length in meters for Reynolds number.
    
    Attributes:
        results (dict): Stores computed atmospheric properties per altitude.
    """
    
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
            self.results[h] = {}
            self.ISA_properties(h, self.results[h])
            self.mach_number(velocity, self.results[h])
            self.reynolds_number(velocity, length, self.results[h])
    
    def ISA_properties(self, altitude, result):
        """
        Computes temperature, pressure, density, and gravity at a given altitude.
        
        Args:
            altitude (float): Altitude in meters.
            result (dict): Dictionary to store results.
        
        Raises:
            ValueError: If altitude is outside the supported range [0, 84852].
        """
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
        """
        Calculates T and P for a gradient (non-isothermal) layer.
        
        Args:
            T (float): Temperature at base of layer.
            P (float): Pressure at base of layer.
            altitude (float): Target altitude.
            h_0 (float): Base altitude of the layer.
            alpha (float): Lapse rate.
        
        Returns:
            Tuple[float, float]: Temperature and Pressure at target altitude.
        """
        T = T + alpha * (altitude - h_0)
        P *= (T / (T - (alpha * (altitude - h_0)))) ** (-self.g0 / (self.R * alpha))
        return T, P
    
    def pause_layer_calculation(self, T, P, altitude, h_0):
        """
        Calculates pressure for an isothermal layer.
        
        Args:
            T (float): Constant temperature of layer.
            P (float): Pressure at base of layer.
            altitude (float): Target altitude.
            h_0 (float): Base altitude of the layer.
        
        Returns:
            float: Pressure at target altitude.
        """
        P *= math.exp((-self.g0 / (self.R * T)) * (altitude - h_0))
        return P
    
    def dynamic_viscosity(self, T):
        """
        Computes dynamic viscosity using Sutherland's formula.
        
        Args:
            T (float): Temperature in Kelvin.
        
        Returns:
            float: Dynamic viscosity [kg/m/s].
        """
        return self.mu_0 * (((T / self.T0)**1.5) * ((self.T0 + self.C) / (T + self.C)))
    
    def mach_number(self, V, result):
        """
        Computes Mach number and speed of sound.
        
        Args:
            V (float): Velocity in m/s.
        
        Returns:
            result (dict): Dictionary to update with Mach and speed of sound.
        """
        a = math.sqrt(self.gamma * self.R * result["Temperature [K]"])
        result["Speed of Sound [m/s]"] = a
        result["Mach Number"] = V / a
    
    def reynolds_number(self, V, L, result):
        """
        Computes Reynolds number.
        
        Args:
            V (float): Velocity in m/s.
            L (float): Characteristic length in meters.
        
        Returns:
            result (dict): Dictionary to update with Reynolds number.
        """
        mu = self.dynamic_viscosity(result["Temperature [K]"])
        result["Reynolds Number"] = (result["Density [kg/m3]"] * V * L) / mu
    
    def __str__(self):
        """
        Creates a formatted string of the ISA results for all altitudes.
        
        Returns:
            str: Formatted output of results.
        """
        output = ["\n"]
        output.append("ISA Calculator Results:")
        for altitude in sorted(self.results):
            res = self.results[altitude]
            output.append(f"Altitude: {altitude:.0f} m")
            max_key_length = max(len(key) for key in res)
            header = f"{'Parameter'.ljust(max_key_length)} | Value"
            output.append(header)
            output.append("-" * len(header))
            for key, value in res.items():
                try:
                    output.append(f"{key.ljust(max_key_length)} | {value:>13,.5f}")
                except (ValueError, TypeError):
                    output.append(f"{key.ljust(max_key_length)} | {'N/A':>13}")
            output.append("")
        return "\n".join(output)
    
    def plot(self):
        """
        Placeholder for future plotting functionality.
        """
        pass


if __name__ == "__main__":
    isa = ISA_Calculator(50000, 50, 2)
    print(isa)