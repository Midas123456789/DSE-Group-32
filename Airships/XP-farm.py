from Airship import Airship
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import fmin


velocity = 25*1.94384 #m/s to kt



hs = np.linspace(10000, 20000, 100)
powers = np.zeros(len(hs))

volumes = np.zeros(len(hs))

cds= np.zeros(len(hs))
airships = {}


fig, axs = plt.subplots(2, 2, figsize=(10, 8))

for n in range(1,5):
    airships[n] = []
    for i in range(len(hs)):
        h = hs[i]
        airship = Airship(2.55, 2e6, n, velocity, h*3.28084,1000,0.95)
        if n == 1:
            airship.BR = 0.98
        airship.iterate_to_exact()
        powers[i] = airship.power_required
        volumes[i] = airship.volume#*0.0283168
        cds[i] = airship.CD
        airships[n]+= [airship]
    #plot volume and power
    axs[0, 0].plot(hs/1000, powers,label=f'{n} lobes')
    axs[0, 1].plot(hs/1000, volumes,label=f'{n} lobes')
    axs[1, 0].plot(hs/1000, cds,label=f'{n} lobes')

#plt.xlabel('Altitude (km)')
axs[0, 0].set_xlabel('Altitude (km)')
axs[0, 1].set_xlabel('Altitude (km)')
axs[0, 0].set_ylabel('Power (kW)')
axs[0, 1].set_ylabel('Volume ($ft^3$)')

axs[1, 0].set_xlabel('Altitude (km)')
axs[1, 0].set_ylabel('Drag Coefficient')
axs[1, 1].set_xlabel('Altitude (km)')
#plt.title('Volume vs Altitude')
axs[0, 0].legend()
axs[0, 1].legend()
axs[1, 0].legend()
plt.show()
print('i am only doing this to test the code')