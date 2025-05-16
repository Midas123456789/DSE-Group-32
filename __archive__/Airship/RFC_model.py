from Power.solar_power_model import Power
import numpy as np

#Values from sources

class RFC:
     def __init__(self, power_model: Power, energy_density = 2.6, vol_energy_density = 1.2, efficiency=0.65): #specific energy in kWh/kg, volemtric specific energy in kWh/m^3 efficency in %
            self.power_model = power_model
            self.energy_density = energy_density
            self.vol_energy_density = vol_energy_density
            self.efficiency = efficiency
        
     def rfc_mass_kg(self):
        deficit = self.power_model.max_deficit()
        if deficit is not None:
            energy_kWh = abs(deficit) / 3.6e6 / self.efficiency  # J → kWh, account for round-trip losses
            return energy_kWh / self.energy_density
        return 0

     def rfc_volume_m3(self):
        deficit = self.power_model.max_deficit()
        if deficit is not None:
            energy_kWh = abs(deficit) / 3.6e6 / self.efficiency  # J → kWh, account for round-trip losses
            return energy_kWh / self.vol_energy_density
        return 0
     
     def required_electrolysis_power(self, electrolyzer_efficiency=0.7):
       deficit = self.power_model.max_deficit()
       if deficit is None:
            return 0

    # Total energy needed to regenerate hydrogen + oxygen (in J)
    # This includes efficiency loss in electrolysis
       energy_needed = abs(deficit) / electrolyzer_efficiency

    # Count the number of seconds during the day when irradiance > 0
       daylight_seconds = np.sum(self.power_model.irradiance > 0)
       if daylight_seconds == 0:
            return 0

    # Required average power during daylight (in Watts)
       return energy_needed / daylight_seconds


     

if __name__ == "__main__":
    # Replace power_required with real data or a test profile
    time = np.arange(86400)
    power_required = 20 * (1 + np.sin(2 * np.pi * time / 86400))  # Example: sine profile

    power_model = Power(latitude=0, day_of_year=172, power_required=power_required)
    rfc = RFC(power_model)

    print("RFC Mass (kg):", rfc.rfc_mass_kg())
    print("RFC Volume (m³):", rfc.rfc_volume_m3())
    print(f"Required Electrolysis Power (W): {rfc.required_electrolysis_power()}")
     
        
    


    

        