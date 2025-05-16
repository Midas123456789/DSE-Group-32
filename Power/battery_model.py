from Power.solar_power_model import *

class Battery:
    def __init__(self, energy_density, density, capasity):
        self.energy_density = energy_density
        self.density = density
        self.capasity = capasity
    
    def battery_mass(self):
        return self.capasity / self.energy_density
    
    def battery_volume(self):
        return self.battery_mass() / self.density


if __name__ == '__main__':
    Pr = 15681 #W
    power_required = [Pr for i in range(86400)]
    power = Power(latitude=20, day_of_year=150, area=200, power_required=power_required)
    power.plot_power_profiles()
    capasity = -power.max_deficit()
    bat = Battery(576000, 2000, capasity)
    print(bat.battery_mass())