from solar_power_model import Power
import numpy as np

#Values from sources

class RFC:
     def __init__(self, power_model: Power, energy_density = 0.65, vol_energy_density = 1.2, efficiency=0.5, energy_density_hydrogen=33.33): #specific energy in kWh/kg, volemtric specific energy in kWh/m^3 efficency in %
            self.power_model = power_model
            self.energy_density = energy_density
            self.vol_energy_density = vol_energy_density
            self.efficiency = efficiency
            self.energy_density_hydrogen = energy_density_hydrogen

     def total_deficit_joules(self):
         """
         Calculate the total energy deficit (in joules) over the entire day.
         """
         net_power = self.power_model.net_power()  # Instantaneous power in W
    
         # Sum the negative power values (deficit) over time in joules (W * seconds)
         deficit_joules = np.sum(np.where(net_power < 0, -net_power, 0))  # Negative power indicates deficit

         return deficit_joules

     def energy_deficit_kWh(self):
   
         deficit_joules = self.total_deficit_joules()  # Total deficit in joules
         deficit_kWh = deficit_joules / 3.6e6  # Convert from Joules to kWh (1 kWh = 3600000 J)
         return deficit_kWh

        
     def rfc_mass_kg(self):
        deficit = self.energy_deficit_kWh()
        energy_kWh = deficit / self.efficiency  
        return energy_kWh / self.energy_density
    

     def rfc_volume_m3(self):
        deficit = self.energy_deficit_kWh()
        energy_kWh = deficit / self.efficiency
        return energy_kWh / self.vol_energy_density
        

     def hydrogen_used(self):
         """
         Calculate the amount of hydrogen used in kg.
         """
         # Get the total energy deficit in kWh
         energy_kWh = self.energy_deficit_kWh() / 0.6
         # Calculate the mass of hydrogen
         mass_hydrogen = energy_kWh / self.energy_density_hydrogen 
         return mass_hydrogen
     
     def water_mass(self):
         mass_water = self.hydrogen_used() * 9
         return mass_water
     
     def required_electrolysis_power(self):
         mass_hydrogen = self.hydrogen_used()
         daylight_seconds = np.sum(self.power_model.irradiance > 0)
         daylight_hours = daylight_seconds / 3600


         energy_required_electrolysis_kWh = mass_hydrogen * 53.4 #kWh
         energy_required_electrolysis_W = energy_required_electrolysis_kWh / daylight_hours * 1000 # Convert kWh to W and divide by seconds in a day

         return energy_required_electrolysis_W
     
     def compute_solar_array_area(self, base_power=100000, efficiency=0.30, solar_constant=1300):
        elevation = self.power_model._solar_elevation()
        raw_irradiance = np.sin(elevation) * solar_constant
        raw_irradiance[elevation <= 0] = 0

        daylight_irradiance = raw_irradiance[raw_irradiance > 0]
        if len(daylight_irradiance) == 0:
            raise ValueError("No daylight irradiance available for solar sizing.")

        avg_irradiance = np.mean(daylight_irradiance)
        electrolysis_power = self.required_electrolysis_power()
        total_power = base_power + electrolysis_power

        area = total_power / (efficiency * avg_irradiance)
        return area
     
     def solar_panel_mass(self):
        """
        Calculate the mass of the solar panels.
        """
        # Assuming a specific mass of 10 kg/m² for solar panels
        specific_mass = 0.8
        area = self.power_model.area
        return area * specific_mass

        

if __name__ == "__main__":
    # Replace power_required with real data or a test profile
    time = np.arange(86400)

    power_required = [100000 for i in range(86400)]

    power_model = Power(latitude=40, day_of_year=1, power_required=power_required, area=30000)
    rfc = RFC(power_model)

    area = rfc.compute_solar_array_area()
    # Step 4: Rebuild the model with the correct profile
    power_model = Power(latitude=40, day_of_year=1, power_required=power_required, area = area)
    rfc = RFC(power_model)

    print("RFC Mass (kg):", rfc.rfc_mass_kg())
    print("RFC Volume (m³):", rfc.rfc_volume_m3())
    print(f"Required Electrolysis Power (W): {rfc.required_electrolysis_power()}")
    print(f"Hydrogen used (kg): {rfc.hydrogen_used()}")
    print(f"Water mass (kg): {rfc.water_mass()}")
    print("areas", area)
    power = power_model.power_generated()
    average_power_W = np.mean(power)
    print(f"Average Power Generated (W): {average_power_W}")
    print(power_model.max_deficit())
    print(rfc.energy_deficit_kWh())
    print(rfc.solar_panel_mass())
    print(rfc.power_model.area)
        
    


    

        