import numpy as np
import matplotlib.pyplot as plt

c = 3e8
k = 1.38e-23  # Boltzmann constant in J/K
T = 290  # Noise temperature in Kelvin
'''
def calculate_link_budget(frequency_mhz, distance_km, tx_power_dbm, tx_gain_db, rx_gain_db, bandwidth_mhz, noise_figure_db):
    # Constants
      # Speed of light in m/s
    

    # Convert frequency to Hz and distance to meters
    frequency_hz = frequency_mhz * 1e6
    distance_m = distance_km * 1e3

    # Free-space path loss (FSPL) in dB
    fspl_db = 20 * math.log10(distance_m) + 20 * math.log10(frequency_hz) - 147.55

    # Noise power in dBm
    noise_power_dbm = 10 * math.log10(k * T * bandwidth_mhz * 1e6) + 30

    # Received power in dBm
    rx_power_dbm = tx_power_dbm + tx_gain_db + rx_gain_db - fspl_db

    # Signal-to-noise ratio (SNR) in dB
    snr_db = rx_power_dbm - noise_power_dbm - noise_figure_db

    return {
        "FSPL (dB)": fspl_db,
        "Noise Power (dBm)": noise_power_dbm,
        "Received Power (dBm)": rx_power_dbm,
        "SNR (dB)": snr_db
    }

#for the haps to cell
# Example parameters for HAPS 5G link
frequency_mhz = 1.6  # 1.6GHz
distance_km = 20  # 20 km
tx_power_dbm = 30  # 1 W
tx_gain_db = 15  # 15 dBi
rx_gain_db = 0  # 15 dBi
bandwidth_mhz = 20  # 20 MHz
noise_figure_db = 5  # 5 dB

# Calculate link budget
link_budget = calculate_link_budget(frequency_mhz, distance_km, tx_power_dbm, tx_gain_db, rx_gain_db, bandwidth_mhz, noise_figure_db)

# Print results
print("Link Budget Results:")
for key, value in link_budget.items():
    print(f"{key}: {value:.2f} dB")

'''
#minimize power for given distance


#maximize ditance for max Power


#maximize snr (capacity) for given distance and X power


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
    #noise += -174+ np.log10(bandwidth)+noise_figure_db # general noise
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
def distance_to_radius(distance,height):
    # Convert distance in km to radius in km
    return distance**2 - height**2

#power_vs_distance(1000, 1.6e9, 0.02,acceptable_snr_db, bandwidth=20e6, noise_figure_db=5) # Example usage

#plot power vs distance
distances = np.linspace(1, 100, 100)  # Distance from 1 km to 100 km
acceptable_snr_db = 10*np.log10(bitrate_to_snr(100e6,20e6))  # Acceptable SNR in dB
print(f'acceptable_snr_db: {acceptable_snr_db} dB')
antenna_diameter = 3  # Antenna diameter in meters
distances*= 1000
powers= [power_vs_distance(distance, 1.6e9, antenna_diameter, acceptable_snr_db, bandwidth=20e6, noise_figure_db=5) for distance in distances]
plt.plot(distances/1000, powers)
plt.xlabel('Distance (km)')
plt.ylabel('Transmitter Power (Watts)')
plt.title('Transmitter Power vs Distance')
plt.grid()
plt.show()