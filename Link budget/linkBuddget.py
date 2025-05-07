import numpy as np
import matplotlib.pyplot as plt

c = 3e8
k = 1.38e-23  # Boltzmann constant in J/K
T = 290  # Noise temperature in Kelvin
height = 20e3  # Height of HAPS in meters. the idea is that this is an input parameter

#minimize power for given distance?


#maximize distance for a set max Power?
#in a way done



#maximize snr (capacity) for given distance and X power?


#graph distancce vs power for  given snr and bandwith
def power_vs_distance(distance,frequency, diameter_antenna,snr_db,bandwidth=20e6, noise_figure_db=5):
    noise_figure = 10**(noise_figure_db/10) # Convert dB to linear scale
    l_fs_db = 10 * np.log10((4 * np.pi * distance * frequency / c)**2) # Free space loss in dB
    l_misc = 2.6 # Miscellaneous losses (e.g., connectors, cables, atmosphere etc.) rnadom web shit
    rx_gain_db = 0 # Receiver gain in dB (assumed to be 0 for simplicity)
   
    tx_gain_db = 10 * np.log10(4*np.pi*0.5*np.pi*(diameter_antenna/2)**2 *(frequency/c)**2)
    #tx_gain_db = 10 * np.log10(4*np.pi*64*0.5)
   
    P_rx_without_tx=tx_gain_db+rx_gain_db - l_fs_db-l_misc
    noise = 10* np.log10(k*(noise_figure-1)*T*bandwidth) # Noise power in dBm
    #noise += -174+ np.log10(bandwidth)+noise_figure_db # general noise in the region and from some paper no clue if correct, though something should be included
    Ptx=noise + snr_db - P_rx_without_tx

    #convert dbm to linear scale
    Ptx = 10**((Ptx)/10) # Convert dBm to linear scale
    #print(f'Ptx: {Ptx} dBm \ntx_gain_db: {tx_gain_db} dB \nrx_gain_db: {rx_gain_db} dB \nFree space loss: {l_fs_db} dB \nMiscellaneous losses: {l_misc} dB \nNoise figure: {noise_figure_db} dB \nNoise power: {noise} dBm \nSNR: {snr_db} dB \nDistance: {distance} m \nFrequency: {frequency} Hz \nAntenna diameter: {diameter_antenna} m')
    

    return Ptx

def bitrate_to_snr(bitrate, bandwidth):
    # Shannon's capacity formula: C = B * log2(1 + SNR)
    # Rearranging gives: SNR = 2^(C/B) - 1
    snr = 2**(bitrate / bandwidth) - 1
    return snr
def distance_to_radius(distance):
    # Convert distance in km to radius in km
    return distance**2 - height**2


total_haps_bitrate = 100e9  # Total bitrate in bps
n_beams = 64  # Number of beams
P_max = 5e3





bandwidth = 100e6  # Bandwidth in Hz; from wikipedia (however not non-teresterial stuff)

#assume worth case and 1/2 of the total bitrate is used in total (between all beams)
acceptable_snr_db = 10*np.log10(bitrate_to_snr(0.5*total_haps_bitrate/n_beams,bandwidth))  # Acceptable SNR in dB
antenna_diameter = 3  # Antenna diameter in meters

#find at which distance the power is equal to P_max/n_beams using root finding
from scipy.optimize import fsolve
def equation_to_solve(distance):
    return power_vs_distance(distance, 3.7e9, antenna_diameter, acceptable_snr_db, bandwidth=bandwidth, noise_figure_db=5) - P_max/n_beams
print(f' the solution is {fsolve(equation_to_solve, 1000)} m') #initial guess of 1000m




# Plotting the graph of power vs distance

#distances = np.linspace(1, 100, 100)  # Distance from 1 km to 100 km
#distances*= 1000
#powers= [power_vs_distance(distance, 3.7e9, antenna_diameter, acceptable_snr_db, bandwidth=20e6, noise_figure_db=5) for distance in distances]
'''
plt.plot(distances/1000, powers)
#do a line at 10kw
plt.axhline(y=P_max/n_beams, color='r', linestyle='--', label='Max Power (10 kW)')
plt.xlabel('Distance (km)')
plt.ylabel('Transmitter Power (Watts)')
plt.title('Transmitter Power vs Distance')
plt.grid()
plt.show()
'''