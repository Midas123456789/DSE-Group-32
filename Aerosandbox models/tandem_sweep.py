from full_model_hydrogen_tandem import get_LH_tandem
from full_model_hydrogen import get_LH_conventional
import matplotlib.pyplot as plt
import numpy as np

W_pl = 1000  # payload weight [N]
P_pl = 10e3  # payload power [W]

def sweep_vs_altitude():
    altitudes = np.linspace(12000, 18000, 6)

    total_masses_conv, total_masses_tand = [], []
    lh_masses_conv, lh_masses_tand = [], []
    power_required_conv, power_required_tand = [], []
    flight_speeds_conv, flight_speeds_tand = [], []

    for altitude in altitudes:
        print(f"Evaluating altitude: {altitude}")
        results_conv = get_LH_conventional(altitude=altitude, P_payload=P_pl, W_payload=W_pl)
        results_tand = get_LH_tandem(altitude=altitude, P_payload=P_pl, W_payload=W_pl)

        total_masses_conv.append(results_conv['W_total'] / 9.81)
        lh_masses_conv.append(results_conv['M_LH'])
        power_required_conv.append(results_conv['P_required'])
        flight_speeds_conv.append(results_conv['V'])

        total_masses_tand.append(results_tand['W_total'] / 9.81)
        lh_masses_tand.append(results_tand['M_LH'])
        power_required_tand.append(results_tand['P_required'])
        flight_speeds_tand.append(results_tand['V'])

    fig, axs = plt.subplots(2, 2, figsize=(12, 10))

    axs[0, 0].plot(altitudes, total_masses_conv, marker='o', label="Conventional")
    axs[0, 0].plot(altitudes, total_masses_tand, marker='s', label="Tandem")
    axs[0, 0].set_title("Total Mass vs Altitude")
    axs[0, 0].set_xlabel("Altitude [m]")
    axs[0, 0].set_ylabel("Total Mass [kg]")
    axs[0, 0].set_ylim(500, 800)
    axs[0, 0].grid(True)
    axs[0, 0].legend()

    axs[0, 1].plot(altitudes, lh_masses_conv, marker='o', label="Conventional", color='green')
    axs[0, 1].plot(altitudes, lh_masses_tand, marker='s', label="Tandem", color='darkgreen')
    axs[0, 1].set_title("LH₂ Mass vs Altitude")
    axs[0, 1].set_xlabel("Altitude [m]")
    axs[0, 1].set_ylabel("LH₂ Mass [kg]")
    axs[0, 1].set_ylim(150, 300)
    axs[0, 1].grid(True)
    axs[0, 1].legend()

    axs[1, 0].plot(altitudes, power_required_conv, marker='o', label="Conventional", color='red')
    axs[1, 0].plot(altitudes, power_required_tand, marker='s', label="Tandem", color='darkred')
    axs[1, 0].set_title("Power Required vs Altitude")
    axs[1, 0].set_xlabel("Altitude [m]")
    axs[1, 0].set_ylabel("Power Required [W]")
    axs[1, 0].set_ylim(2000, 7000)
    axs[1, 0].grid(True)
    axs[1, 0].legend()

    axs[1, 1].plot(altitudes, flight_speeds_conv, marker='o', label="Conventional", color='purple')
    axs[1, 1].plot(altitudes, flight_speeds_tand, marker='s', label="Tandem", color='indigo')
    axs[1, 1].set_title("Flight Speed vs Altitude")
    axs[1, 1].set_xlabel("Altitude [m]")
    axs[1, 1].set_ylabel("Speed [m/s]")
    axs[1, 1].set_ylim(10, 40)
    axs[1, 1].grid(True)
    axs[1, 1].legend()

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    sweep_vs_altitude()
