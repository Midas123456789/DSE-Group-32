from Airship import Airship
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import fmin

max_br_hybrid = 0.85
max_br_conv = 0.95
min_br=0.7
FR_range = [1.2, 15]
velocity = 20*1.94384 #m/s to kt


def calculate_power(airship, vec):
    BR = vec[0]
    FR = vec[1]
    if BR<min_br or BR>max_br_hybrid:
        if airship.n_lobes == 1 and BR > max_br_conv:
            return 100000000000000
        elif airship.n_lobes != 1:
            return 100000000000000
    if FR<FR_range[0] or FR>FR_range[1]:
        return 100000000000000
    airship.BR = BR
    airship.FR = FR
    airship.iterate_to_exact()
    #return (airship.CD0 + airship.K * airship.CL ** 2)*airship.q*airship.reference_volume *airship.velocity
    return airship.power_required


hs = np.linspace(10000, 20000, 100)
powers = np.zeros(len(hs))
brs = np.zeros(len(hs))
frs = np.zeros(len(hs))
volumes = np.zeros(len(hs))

airships = {}

fig, axs = plt.subplots(2, 2, figsize=(10, 8))

for n in range(1,5):
    airships[n] = []
    for i in range(len(hs)):
        h = hs[i]
        airship = Airship(3, 9e6, n, velocity, h*3.28084,1000,0.7)
        p = lambda vec: calculate_power(airship, vec)
        #airship.iterate_to_exact()
        #print(airship)

        fmin(p, [0.8,3], disp=False) #[BR,FR]
        #fmin(v, 0.8, disp=False)
        #powers[i] = (airship.CD0 + airship.K * airship.CL ** 2)*airship.q*airship.reference_volume *airship.velocity
        powers[i] = airship.power_required
        brs[i] = airship.BR
        volumes[i] = airship.volume#*0.0283168
        frs[i] = airship.FR
        airships[n]+= [airship]
    #BRs = np.linspace(0.8, 0.95, 100)
    #airship = Airship(3, 2e8, 3, 84.5, 60e3, 1000, 0.7)
    #make 3 graphs, one with volume, one with power, one with BR

    axs[0, 0].plot(hs/1000, powers,label=f'{n} lobes')
    axs[0, 1].plot(hs/1000, brs,label=f'{n} lobes')
    axs[1, 0].plot(hs/1000, volumes,label=f'{n} lobes')
    axs[1, 1].plot(hs/1000, frs,label=f'{n} lobes')


axs[0, 0].set_xlabel('Altitude (km)')
axs[0, 0].set_ylabel('Power (kW)')
axs[0, 0].set_title('Power vs Altitude')
axs[0, 0].grid()
axs[0, 0].legend()

axs[0, 1].set_xlabel('Altitude (km)')
axs[0, 1].set_ylabel('BR')
axs[0, 1].set_title('BR vs Altitude')
axs[0, 1].grid()
axs[0, 1].legend()

axs[1, 0].set_xlabel('Altitude (km)')
axs[1, 0].set_ylabel('Volume ($ft^3$)')
axs[1, 0].set_title('Volume vs Altitude')
axs[1, 0].grid()
axs[1, 0].legend()

axs[1, 1].set_xlabel('Altitude (km)')
axs[1, 1].set_ylabel('FR')
axs[1, 1].set_title('FR vs Altitude')
axs[1, 1].grid()
axs[1, 1].legend()
plt.show()

#clean up the figure
plt.figure(figsize=(10, 5))
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

#plot the difference between wg1 and wg2 for each airship as a function of altitude


diffs =[ [airship.wg - airship.W_g2 for airship in airships[n]] for n in range(1,5)]
for i in range(4):
    plt.plot(hs/1000, diffs[i], label=f'{i+1} lobes')
plt.xlabel('Altitude (km)')
plt.ylabel('Difference in Weight (N)')
plt.title('Difference in Weight vs Altitude')
plt.grid()
plt.legend()
plt.show()

#for airships with 1 lobe, pick airhips with difference in weight more than 0.01
problematic_airships = []
for i in range(len(airships[1])):
    if abs(airships[1][i].wg - airships[1][i].W_g2) > 0.01:
        problematic_airships+=[airships[1][i]]

#pick the first 5 airships
problematic_airships = problematic_airships[:5]
#plot the difference in weight for these airships as a function of volume
plt.figure(figsize=(10, 5))
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)



for i in range(len(problematic_airships)):
    airship = problematic_airships[i]
    vols= np.linspace(airship.volume*0.5, 2e7,100)
    diffs = np.zeros(len(vols))
    for i in range(len(vols)):
        
        diffs[i] = airship.iterator([vols[i]])
    plt.plot(vols, diffs, label=f'altitude {airship.altitude/(3.28084*1000):.2f} km')
plt.xlabel('Volume ($ft^3$)')
plt.ylabel('Difference in Weight (lb)')
plt.title('Difference in Weight vs Volume')
plt.grid()
plt.legend()
plt.ylim(-3e5, 3e5)
plt.show()
print('idk')
