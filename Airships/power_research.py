#see what max weight of rfc subsystem  we can have for the model to still converge

#W=>P=>

#18km has the least wind
from Airship import Airship
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import fmin
from RFC_model import RFC
from solar_power_model import Power
from scipy.stats import linregress

#fix altitude, velocity, and n_lobes
velocity = 25*1.94384 #m/s to kt
h = 20*3.28084 #ft

def make_rfc_to_line():
    #see how the RFC behaves with different powers
    powers = np.linspace(1e5, 1e8, 100)
    rfcs = []

    for p in powers:
        power_model = Power(latitude=55, day_of_year=150, power_required=p, area= 1000)
        rfc = RFC(power_model)

        compt_area = rfc.compute_solar_array_area()
        # Step 4: Rebuild the model with the correct profile
        power_model = Power(latitude=55, day_of_year=150, power_required=p, area = compt_area)
        rfcs += [RFC(power_model)]


    #do linear fit for the mass of the RFC
    # Calculate total mass for each RFC
    masses = [rfc.rfc_mass_kg() + rfc.solar_panel_mass() for rfc in rfcs]

    # Perform linear regression
    slope, intercept, r_value, p_value, std_err = linregress(powers, masses)
    print(f"Slope: {slope}, Intercept: {intercept}, R-squared: {r_value**2}")
    # Plot the data and the linear fit
    plt.plot(powers, masses, label='RFC Mass')
    plt.plot(powers, slope * powers + intercept, 'r--', label='Linear Fit')
    plt.legend()
    #plt.plot(powers, 4.14e-2 *powers, label='RFC Mass')
    plt.xlabel('Power (W)')
    plt.ylabel('Mass (kg)')
    plt.title('RFC Mass vs Power required')
    plt.grid()
    plt.show()

#find the max slope for the airship that still ensures convergence per altitude

airship = Airship(3, 3e6, 1, velocity, h, 250**2.20462, 0.95)
airship.iterate_to_exact()
print(airship)
xx= airship.iterator([airship.volume])
print(airship.power_required)