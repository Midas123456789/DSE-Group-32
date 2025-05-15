from Airship import Airship
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import fmin


velocity = 25*1.94384 #m/s to kt



hs = np.linspace(1000, 20000, 100)
powers = np.zeros(len(hs))

volumes = np.zeros(len(hs))

cds= np.zeros(len(hs))
airships = {}


fig, axs = plt.subplots(3, 3, figsize=(10, 8))

for n in range(1,5):
    airships[n] = []
    for i in range(len(hs)):
        h = hs[i]
        airship = Airship(3, 2e7, n, velocity, h*3.28084,225,0.8)
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
    axs[1, 1].plot(hs/1000, [airship.iterator([airship.volume]) for airship in airships[n]],label=f'{n} lobes')

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
#plt.clf()

#grpaph the weight difference as a finction of volume from 10^5 to 10^7 for all the airships
vols = np.linspace(1e5, 1e8, 100)
for n in range(1,5):
    airship =  airships[n][0]
    plt.plot(vols, [airship.iterator([vol]) for vol in vols] , label=f'{n} lobes')
    plt.xlabel('Volume ($ft^3$)')
    plt.ylabel('Weight diff (lb)')
    plt.ylim(-1e5, 100)
plt.legend()
plt.grid()
plt.show()