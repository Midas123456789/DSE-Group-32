import numpy as np
import matplotlib.pyplot as plt

class PropPowerAirship:
    def __init__(self, airship_CD0, drag_factor_K, airship_CLaero, airship_velocity, rho, airship_refrencearea, prop_efficiency=0.8):
    
        self.airship_CD0 = airship_CD0
        self.drag_factor_K = drag_factor_K
        self.airship_CLaero = airship_CLaero
        self.airship_velocity = airship_velocity
        self.rho = rho
        self.airship_refrencearea = airship_refrencearea
        self.prop_efficiency = prop_efficiency
        self.payload = 2.75       #KW

    
    def prop_power(self):
        prop_power_req = (self.airship_velocity*((self.airship_CD0 + self.drag_factor_K * self.airship_CLaero**2) * 0.5 * self.rho * self.airship_velocity**2 * self.airship_refrencearea)) / (550*(self.prop_efficiency)) #ft-lb/s
        self.prop_power_req = prop_power_req * 0.745699872 #* (1.35582/1000)  # Convert from ft-lb/s to kWatts
        return self.prop_power_req

    def total_power(self):

        self.total_power = self.prop_power() + self.payload
        return self.total_power

"""
if __name__ == "__main__":
    # Constants
    airship_CD0 = 0.025 
    drag_factor_K = 0.05
    airship_CLaero = 0.2
    rho = 0.00017266  # slugs/ft^3
    airship_refrencearea = 10562   # ft^2
    prop_efficiency = 0.63

    # Airship velocity range
    airship_velocity = np.linspace(10, 150, 100)  # ft/s

    # Create an instance of the PropPowerAirship class
    prop_power_airship = PropPowerAirship(airship_CD0, drag_factor_K, airship_CLaero, airship_velocity, rho, airship_refrencearea, prop_efficiency)

    # Calculate power required
    power_required = prop_power_airship.calculate_power()

    # Plotting the results
    plt.plot(airship_velocity*0.3048, power_required)
    plt.title('Power Required vs Airship Velocity')
    plt.xlabel('Airship Velocity (ms/s)')
    plt.ylabel('Power Required (kW)')
    plt.grid()
    plt.show()
"""


