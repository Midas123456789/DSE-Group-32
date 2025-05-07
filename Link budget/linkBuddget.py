import math

def calculate_link_budget(frequency_mhz, distance_km, tx_power_dbm, tx_gain_db, rx_gain_db, bandwidth_mhz, noise_figure_db):
    # Constants
    c = 3e8  # Speed of light in m/s
    k = 1.38e-23  # Boltzmann constant in J/K
    T = 290  # Noise temperature in Kelvin

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


#minimize power for given distance


#maximize ditance for max Power


#maximize snr (capacity) for given distance and X power


#graph distancce vs power for  given snr and bandwith
def power_vs_distance(frequency_mhz, tx_power_dbm, tx_gain_db, rx_gain_db, bandwidth_mhz, noise_figure_db):
    distance
    powers = []

    for distance_km in distances:
        link_budget = calculate_link_budget(frequency_mhz, distance_km, tx_power_dbm, tx_gain_db, rx_gain_db, bandwidth_mhz, noise_figure_db)
        powers.append(link_budget["Received Power (dBm)"])

    return powers