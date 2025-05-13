import aerosandbox as asb
import aerosandbox.numpy as np
from mass_wing import Mass_wing
import csv
import matplotlib.pyplot as plt

solar_irradiance_values = np.linspace(300, 1000, 20)
results = []

wing_airfoil = asb.Airfoil("sd7037")
N_cords = 4
g = 9.81
altitude = 16000
atm = asb.Atmosphere(altitude=altitude)

solar_efficiency = 0.20
motor_eff = 0.8
prop_eff = 0.8
η_prop = motor_eff * prop_eff
Tail_drag_factor = 1.05

day_length = 0.5
second_in_day = 86400
Energy_density_bat = 0.9e6  # J/kg

skin_density = 0.25  # kg/m²
W_payload = 580  # N
P_payload = 5000

for irradiance in solar_irradiance_values:
    print('==========================', irradiance, '=======================')
    try:
        opti = asb.Opti()
        cords = opti.variable(init_guess=3 * np.ones(N_cords), n_vars=N_cords)
        b = opti.variable(init_guess=60, upper_bound=200, lower_bound=0)
        V = opti.variable(init_guess=30, upper_bound=200, lower_bound=1)
        alpha = opti.variable(init_guess=5, lower_bound=0, upper_bound=30)

        y_sections = np.linspace(0, b / 2, N_cords)
        wing = Mass_wing(
            symmetric=True,
            xsecs=[
                asb.WingXSec(
                    xyz_le=[-0.25 * cords[i], y_sections[i], 0],
                    chord=cords[i],
                    airfoil=wing_airfoil
                )
                for i in range(N_cords)
            ]
        ).translate([4, 0, 0])
        airplane = asb.Airplane(wings=[wing])

        opti.subject_to([
            cords > 0,
            np.diff(cords) <= 0,
            wing.area() < 2000
        ])

        op = asb.OperatingPoint(velocity=V, alpha=alpha, atmosphere=atm)
        aero = asb.AeroBuildup(airplane=airplane, op_point=op).run()

        L = 0.5 * atm.density() * V**2 * aero['CL'] * wing.area()
        D = 0.5 * atm.density() * V**2 * aero['CD'] * wing.area() * Tail_drag_factor

        P_required = V * D
        E_generated = wing.area() * irradiance * solar_efficiency * day_length * second_in_day
        E_payload = P_payload * second_in_day
        E_available_for_prop = (E_generated - E_payload) * η_prop
        E_required = P_required * second_in_day

        opti.subject_to(E_available_for_prop >= E_required)

        Bat_mass = E_required / Energy_density_bat
        W_battery = Bat_mass * g
        W_skin = 2 * wing.area() * skin_density * g
        Skin_mass = 2 * wing.area() * skin_density
        spar_mass, t_list = wing.spar_mass(L, b)
        W_spar = spar_mass * g
        W_wing = W_spar + W_skin
        W_rest = W_wing
        W_total = W_payload + W_wing + W_battery + 0*W_rest
        opti.subject_to(L >= W_total)

        opti.minimize(b)
        sol = opti.solve()

        results.append({
            "Irradiance": irradiance,
            "Span": sol.value(b),
            "Wing Area": sol.value(wing.area()),
            "Battery Mass": sol.value(Bat_mass),
            "Skin Mass": sol.value(Skin_mass),
            "Spar Mass": sol.value(spar_mass),
            "Structural Mass": sol.value(Skin_mass + spar_mass),
            "Total Weight": sol.value(W_total),
            "Velocity": sol.value(V),
        })

    except Exception as e:
        print(f"Optimization failed for irradiance {irradiance} W/m²: {e}")

# === Postprocess === #
for row in results:
    row["Wing Loading"] = row["Total Weight"] / row["Wing Area"]

# === Save to CSV === #
with open("solar_sweep_results.csv", "w", newline="") as csvfile:
    fieldnames = ["Irradiance", "Span", "Wing Area", "Battery Mass", "Skin Mass", "Spar Mass", "Structural Mass", "Total Weight", "Velocity", "Wing Loading"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for row in results:
        writer.writerow(row)

print("Results saved to solar_sweep_results.csv")

# === Extract Data === #
irr_vals = [r["Irradiance"] for r in results]
span_vals = [r["Span"] for r in results]
area_vals = [r["Wing Area"] for r in results]
bat_mass_vals = [r["Battery Mass"] for r in results]
skin_mass_vals = [r["Skin Mass"] for r in results]
spar_mass_vals = [r["Spar Mass"] for r in results]
struct_mass_vals = [r["Structural Mass"] for r in results]
weight_vals = [r["Total Weight"] for r in results]
velocity_vals = [r["Velocity"] for r in results]
wing_loading_vals = [r["Wing Loading"] for r in results]

# === Plot Design Results === #
plt.figure(figsize=(10, 14))

def make_plot(index, y, ylabel, title):
    ax = plt.subplot(3, 2, index)
    ax.plot(irr_vals, y, marker='o')
    ax.set_xlabel("Solar Irradiance (W/m²)")
    ax.set_ylabel(ylabel)
    ax.grid(True)
    ax.text(irr_vals[0] + 20, max(y) - 0.05 * (max(y) - min(y)), title,
            fontsize=10, fontweight='bold',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

make_plot(1, span_vals, "Wingspan (m)", "Wingspan vs Irradiance")
make_plot(2, area_vals, "Wing Area (m²)", "Wing Area vs Irradiance")
make_plot(3, bat_mass_vals, "Battery Mass (kg)", "Battery Mass vs Irradiance")
make_plot(4, weight_vals, "Total Weight (N)", "Total Weight vs Irradiance")
make_plot(5, velocity_vals, "Velocity (m/s)", "Velocity vs Irradiance")
make_plot(6, wing_loading_vals, "Wing Loading (N/m²)", "Wing Loading vs Irradiance")

plt.tight_layout()
plt.show()

# === Plot Subsystem Masses === #
plt.figure(figsize=(10, 6))

plt.plot(irr_vals, bat_mass_vals, marker='o', label="Battery Mass")
plt.plot(irr_vals, skin_mass_vals, marker='s', label="Skin Mass")
plt.plot(irr_vals, spar_mass_vals, marker='^', label="Spar Mass")
plt.plot(irr_vals, struct_mass_vals, marker='x', label="Structural Mass")

plt.xlabel("Solar Irradiance (W/m²)")
plt.ylabel("Mass (kg)")
plt.title("Subsystem Masses vs Solar Irradiance")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
