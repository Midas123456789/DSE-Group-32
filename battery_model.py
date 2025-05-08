from solar_power_model import *

class Battery:
    def __init__(self, energy_density, density, capasity):
        self.energy_density = energy_density
        self.density = density
        self.capasity = capasity
    
    def battery_mass(self):
        return self.capasity * self.energy_density
    
    def battery_volume(self):
        return self.battery_mass() / self.density
