import numpy as np
import matplotlib.pyplot as plt
class Power:
    def __init__(self, latitude: float, day_of_year: int, power_required=None, max_irradiance: float = 1000, efficiency: float = 0.2):
        self.latitude = np.radians(latitude)
        self.day_of_year = day_of_year
        self.max_irradiance = max_irradiance
        self.efficiency = efficiency
        self.seconds_in_day = 86400
        self.seconds = np.arange(self.seconds_in_day)
        self.irradiance = self._compute_irradiance()
        self.power_required = (
            np.array(power_required) if power_required is not None
            else np.zeros(self.seconds_in_day)
        )

    def _declination(self):
        return np.radians(23.45) * np.sin(np.radians(360 / 365 * (284 + self.day_of_year)))

    def _hour_angle(self):
        solar_time = self.seconds / 3600  # in hours
        ha_deg = (solar_time - 12) * 15
        return np.radians(ha_deg)

    def _solar_elevation(self):
        decl = self._declination()
        ha = self._hour_angle()
        return np.arcsin(
            np.sin(self.latitude) * np.sin(decl) +
            np.cos(self.latitude) * np.cos(decl) * np.cos(ha)
        )

    def _compute_irradiance(self):
        elevation = self._solar_elevation()
        raw_irradiance = self.max_irradiance * np.sin(elevation)
        raw_irradiance[elevation <= 0] = 0
        return raw_irradiance

    def power_generated(self):
        return self.irradiance * self.efficiency

    def net_power(self):
        return self.power_generated() - self.power_required

    def net_energy(self):
        return np.cumsum(self.net_power())

    def max_surplus(self):
        net = self.net_energy()
        return net.max() if net.max() > 0 else None

    def max_deficit(self):
        net = self.net_energy()
        return net.min() if net.min() < 0 else None
    
    def plot_power_profiles(self):
        seconds = self.seconds
        generated = self.power_generated()
        required = self.power_required
        net = self.net_power()
        energy = self.net_energy()

        fig, ax1 = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

        # Plot power values
        ax1[0].plot(seconds, generated, label='Power Generated (W)', color='green')
        ax1[0].plot(seconds, required, label='Power Required (W)', color='red', linestyle='--')
        ax1[0].plot(seconds, net, label='Net Power (W)', color='blue')
        ax1[0].set_ylabel('Power (W)')
        ax1[0].legend()
        ax1[0].grid(True)
        ax1[0].set_title('Power Profiles Throughout the Day')

        # Plot net energy
        ax1[1].plot(seconds, energy, label='Net Energy (J)', color='purple')
        ax1[1].set_ylabel('Energy (J)')
        ax1[1].set_xlabel('Seconds Since Midnight')
        ax1[1].legend()
        ax1[1].grid(True)
        ax1[1].set_title('Cumulative Net Energy')

        plt.tight_layout()
        plt.show()




if __name__ == "__main__":
    seconds_in_day = 86400
    time = np.arange(seconds_in_day)
    amplitude = 20  # Max power required in watts
    power_required = amplitude * (1 + np.sin(2 * np.pi * time / seconds_in_day))  # Sine wave with period of 24 hours

    power = Power(latitude=0, day_of_year=0, power_required=power_required)

    print('Deficit: ', power.max_deficit())
    print('Surplus: ', power.max_surplus())



    power.plot_power_profiles()


