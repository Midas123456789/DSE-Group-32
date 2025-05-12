from airship_customized import Airship
import numpy as np


hs = np.linspace(10000, 20000, 100)
cls=np.zeros(len(hs))
vols = np.zeros(len(hs))
for i in range(len(hs)):
    h = hs[i]
    airship = Airship(3, 2e8, 3, 84.5, h*3.28084,1000,0.7)
    airship.iterate_to_exact()
    print(airship)
    vols[i] = airship.volume
    cls[i] = airship.CL
import matplotlib.pyplot as plt
plt.plot(hs/1000, vols)
plt.xlabel('Altitude (km)')
plt.ylabel('Volume (m^3)')
plt.grid()
plt.show()

plt.plot(hs/1000, cls)
plt.xlabel('Altitude (km)')
plt.ylabel('$C_L$')
plt.grid()
plt.show()