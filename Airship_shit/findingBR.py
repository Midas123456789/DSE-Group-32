from airship_customized import Airship
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import fmin
max_br = 0.98
def calculate_power(airship, BR):
    if BR<0.8 or BR>max_br:
        return 100000000000000
    airship.BR = BR[0]
    airship.iterate_to_exact()
    return (airship.CD0 + airship.K * airship.CL ** 2)*airship.q*airship.reference_volume *airship.velocity
def calculate_volume(airship, BR):
    if BR<0.8 or BR>max_br:
        return 100000000000000
    airship.BR = BR[0]
    airship.iterate_to_exact()
    return airship.volume
hs = np.linspace(10000, 20000, 100)
powers = np.zeros(len(hs))
brs = np.zeros(len(hs))
volumes = np.zeros(len(hs))
for i in range(len(hs)):
    h = hs[i]
    airship = Airship(3, 2e6, 3, 64.5, h*3.28084,1000,0.7)
    p = lambda BR: calculate_power(airship, BR)
    v = lambda BR: calculate_volume(airship, BR)
    #airship.iterate_to_exact()
    #print(airship)
    
    fmin(p, 0.8, disp=False)
    #fmin(v, 0.8, disp=False)
    powers[i] = (airship.CD0 + airship.K * airship.CL ** 2)*airship.q*airship.reference_volume *airship.velocity
    brs[i] = airship.BR
    volumes[i] = airship.volume
#BRs = np.linspace(0.8, 0.95, 100)
#airship = Airship(3, 2e8, 3, 84.5, 60e3, 1000, 0.7)
#make 3 graphs, one with volume, one with power, one with BR
fig, axs = plt.subplots(2, 2, figsize=(10, 8))
axs[0, 0].plot(hs/1000, powers)
axs[0, 0].set_xlabel('Altitude (km)')
axs[0, 0].set_ylabel('Power (W)')
axs[0, 0].set_title('Power vs Altitude')
axs[0, 0].grid()
axs[0, 1].plot(hs/1000, brs)
axs[0, 1].set_xlabel('Altitude (km)')
axs[0, 1].set_ylabel('BR')
axs[0, 1].set_title('BR vs Altitude')
axs[0, 1].grid()
axs[1, 0].plot(hs/1000, volumes)
axs[1, 0].set_xlabel('Altitude (km)')
axs[1, 0].set_ylabel('Volume (m^3)')
axs[1, 0].set_title('Volume vs Altitude')
axs[1, 0].grid()
plt.show()
