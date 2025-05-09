import numpy as np
from scipy.integrate import quad

class Power:
    def __init__(self, I_max=1000, latitude=52.0, declination=23.44):
        self.I_max = I_max  # Maximum solar irradiance (W/m²)
        self.latitude = latitude  # Observer's latitude (degrees)
        self.declination = declination  # Solar declination (degrees)
    
    def solar_angle(self, t):
        """
        Calculate the solar zenith angle based on time `t` (in hours).
        """
        hour_angle = 15 * (t - 12)  # Hour angle: 15° per hour, centered on noon
        angle = 90 - (self.latitude + self.declination * np.cos(np.deg2rad(hour_angle)))
        return np.deg2rad(angle)
    
    def irradiance(self, t):
        """
        Calculate the solar irradiance at time `t` (in hours) using a cosine model.
        """
        angle = self.solar_angle(t)
        return self.I_max * np.cos(angle)
    
    def energy_over_day(self, sunrise_time=6, sunset_time=18):
        """
        Calculate the total solar energy (in kWh) using numerical integration.
        sunrise_time and sunset_time are in hours.
        """
        # Perform the integration of irradiance over the daylight hours
        result, error = quad(self.irradiance, sunrise_time, sunset_time)
        # Convert energy from Watt-seconds to kWh (1 W = 1 J/s, 1 hour = 3600 seconds)
        energy_kWh = result / 1000 / 3600
        return energy_kWh

# Usage example
power_model = Power(I_max=1000, latitude=52.0, declination=23.44)

# Calculate the total energy from 6 AM to 6 PM
total_energy = power_model.energy_over_day(sunrise_time=6, sunset_time=18)
print(f"Total energy received during the day: {total_energy:.2f} kWh")
