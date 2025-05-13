from airship_customized import Airship
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

#airship = Airship(3, 2e6, 3, 84.5, 4000, 40000, 0.7)
#airship.iterate_to_exact()


"""
hs = np.linspace(10000, 20000, 100)
cls=np.zeros(len(hs))
vols = np.zeros(len(hs))
for i in range(len(hs)):
    h = hs[i]
    airship = Airship(3, 2e8, 3, 84.5, h*3.28084,1000,0.7)
    airship.iterate_to_exact()
    print(airship)

    weight[i] = airship.wg
    vols[i] = airship.volume
    cls[i] = airship.CL
    cds[i] = airship.CD


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
"""


#def altitude_graph(self):
alt_step = 1000
graph_df = [0,0,0,0,0,0,0]
altitudes = np.linspace(10000, 20000, alt_step)

for altitude in altitudes:
    airship = Airship(3,2e5, 3, 84.5, altitude * 3.28084, 1000, 0.7)
    airship.iterate_to_exact()
    print(airship)
    row = [altitude, airship.volume, airship.wg, airship.CL, airship.CD, (airship.buoyancy/airship.D), (airship.CL/airship.CD)]

    graph_df = np.vstack((graph_df, row))

graph_df = pd.DataFrame(graph_df, columns=['altitude', 'Volume', 'WeightG', 'CL', 'CD', 'Lb/D', 'CL/CD'])
graph_df = graph_df.iloc[2:]

fig, axes = plt.subplots(2, 3, figsize=(8, 6))
axes = axes.flatten()  # Flatten the 2x2 array of axes for easy indexing

# Create subplots: 2 rows, 2 columns
fig, axs = plt.subplots(2, 3, figsize=(10, 8))

# Plot in each subplot
axs[0, 0].plot(graph_df['Volume'], graph_df['altitude'])
axs[0, 0].set_title('altitude vs Volume')

axs[0, 1].plot(graph_df['WeightG'], graph_df['altitude'])
axs[0, 1].set_title('altitude vs WeightG')

axs[1, 0].plot(graph_df['CL'], graph_df['altitude'])
axs[1, 0].set_title('altitude vs CL')

axs[1, 1].plot(graph_df['CD'], graph_df['altitude'])
axs[1, 1].set_title('altitude vs CD')

axs[0, 2].plot(graph_df['CL/CD'], graph_df['altitude'])
axs[0, 2].set_title('altitude vs CL/CD')

axs[1, 2].plot(graph_df['Lb/D'], graph_df['altitude'])
axs[1, 2].set_title('altitude vs Lb/D')

# Adjust layout
plt.tight_layout()
plt.show()

#return

