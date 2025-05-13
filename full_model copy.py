import aerosandbox as asb
import aerosandbox.numpy as np
from mass_wing import Mass_wing

# === Constants === #
wing_airfoil = asb.Airfoil("sd7037")
N_cords = 4
g = 9.81
altitude = 16000
atm = asb.Atmosphere(altitude=altitude)
opti = asb.Opti()

# === Solar & Power Parameters === #
solar_irradiance = 700         # W/m^2 available at altitude
solar_efficiency = 0.20        # 20% cell efficiency
motor_eff = 0.8
prop_eff = 0.8
η_prop = motor_eff * prop_eff  # Total propulsive efficiency

# === Mass Parameters === #
W_payload = 800  # N

# === Design Variables === #
cords = opti.variable(init_guess=3 * np.ones(N_cords), n_vars=N_cords)
b = opti.variable(init_guess=60, upper_bound=70, lower_bound=0)
V = opti.variable(init_guess=30, upper_bound=200, lower_bound=1)
alpha = opti.variable(init_guess=5, lower_bound=0, upper_bound=30)

# === Geometry Setup === #
y_sections = np.linspace(0, b / 2, N_cords)

# Main wing
wing_main = Mass_wing(
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

# Second wing (e.g. tailplane or canard)
wing_second = Mass_wing(
    symmetric=True,
    xsecs=[
        asb.WingXSec(
            xyz_le=[-0.25 * cords[i], y_sections[i], 0],
            chord=cords[i],
            airfoil=wing_airfoil
        )
        for i in range(N_cords)
    ]
).translate([10, 0, 0])  # Different location along x-axis

# Full airplane with two wings
airplane = asb.Airplane(wings=[wing_main, wing_second])

# === Constraints === #
opti.subject_to([
    cords > 0,
    np.diff(cords) <= 0,  # Taper
    (wing_main.area() + wing_second.area()) < 2000
])

# === Aerodynamics === #
op = asb.OperatingPoint(velocity=V, alpha=alpha, atmosphere=atm)
aero = asb.AeroBuildup(airplane=airplane, op_point=op).run()

# === Lift & Drag Forces === #
total_area = wing_main.area() + wing_second.area()
L = 0.5 * atm.density() * V**2 * aero['CL'] * total_area
D = 0.5 * atm.density() * V**2 * aero['CD'] * total_area

# === Payload Power === #
P_payload = 10000  # Watts, constant power for onboard systems

# === Power Constraint and Battery Sizing === #
day_length = 0.5
second_in_day = 86400
P_required = V * D
E_generated = total_area * solar_irradiance * solar_efficiency * day_length * second_in_day
E_payload = P_payload * second_in_day
E_available_for_prop = (E_generated - E_payload) * η_prop
E_required = P_required * second_in_day

opti.subject_to(E_available_for_prop >= E_required)

Energy_density_bat = 0.72 * 10**6  # J/kg
Bat_mass = E_required / Energy_density_bat
W_battery = Bat_mass * g

# === Skin Weight === #
skin_density = 0.25  # kg/m²
W_skin = 2 * total_area * skin_density * g

# === Spar Weight === #
spar_mass_main, t_list_main = wing_main.spar_mass(L / 2, b)
spar_mass_second, t_list_second = wing_second.spar_mass(L / 2, b)
W_spar = (spar_mass_main + spar_mass_second) * g

# === Total Weight === #
W_total = W_payload + W_spar + W_skin + W_battery
opti.subject_to(L >= W_total)

# === Objective === #
opti.minimize(b)

# === Solve === #
sol = opti.solve()

# === Extract Solved Values === #
cords_sol = sol.value(cords)
b_sol = sol.value(b)
V_sol = sol.value(V)
alpha_sol = sol.value(alpha)
area_main_sol = sol.value(wing_main.area())
area_second_sol = sol.value(wing_second.area())
area_sol = area_main_sol + area_second_sol
CL_sol = sol.value(aero['CL'])
CD_sol = sol.value(aero['CD'])
L_sol = sol.value(L)
D_sol = sol.value(D)
P_required_sol = sol.value(P_required)
W_spar_sol = sol.value(W_spar)
t_list_main_sol = sol(t_list_main)
t_list_second_sol = sol(t_list_second)
battery_mass_sol = sol(Bat_mass)
Total_weight_sol = sol(W_total)

# === Rebuild Solved Wings === #
y_sections_sol = np.linspace(0, b_sol / 2, N_cords)

wing_main_sol = Mass_wing(
    symmetric=True,
    xsecs=[
        asb.WingXSec(
            xyz_le=[-0.25 * cords_sol[i], y_sections_sol[i], 0],
            chord=cords_sol[i],
            airfoil=wing_airfoil
        )
        for i in range(N_cords)
    ]
).translate([4, 0, 0])

wing_second_sol = Mass_wing(
    symmetric=True,
    xsecs=[
        asb.WingXSec(
            xyz_le=[-0.25 * cords_sol[i], y_sections_sol[i], 0],
            chord=cords_sol[i],
            airfoil=wing_airfoil
        )
        for i in range(N_cords)
    ]
).translate([10, 0, 0])

airplane_sol = asb.Airplane(wings=[wing_main_sol, wing_second_sol])
airplane_sol.draw(show=True)

# === Results === #
print("\n" + "="*35)
print("        OPTIMIZATION RESULTS")
print("="*35 + "\n")

print("=== Geometry ===")
print(f"  Span (b):                 {b_sol:.2f} m")
print(f"  Chord distribution (m):   {cords_sol}")
print(f"  Wing area (main):         {area_main_sol:.2f} m²")
print(f"  Wing area (second):       {area_second_sol:.2f} m²")
print(f"  Total wing area:          {area_sol:.2f} m²\n")

print("=== Flight Conditions ===")
print(f"  Flight velocity (V):      {V_sol:.2f} m/s")
print(f"  Angle of attack (α):      {alpha_sol:.2f}°")
print(f"  Altitude:                 {altitude} m")
print(f"  Air density:              {atm.density():.4f} kg/m³\n")

print("=== Aerodynamics ===")
print(f"  Lift coefficient (CL):    {CL_sol:.3f}")
print(f"  Drag coefficient (CD):    {CD_sol:.3f}")
print(f"  Lift (L):                 {L_sol:.2f} N")
print(f"  Drag (D):                 {D_sol:.2f} N\n")

print("=== Energy & Power ===")
print(f"  Payload power:            {P_payload:.2f} W")
print(f"  Power required (P):       {P_required_sol:.2f} W")
print(f"  Energy required:          {sol.value(E_required) / 1e6:.2f} MJ")
print(f"  Energy available (solar): {sol.value(E_available_for_prop) / 1e6:.2f} MJ")
print(f"  Solar efficiency:         {solar_efficiency*100:.1f}%")
print(f"  Propulsion efficiency:    {η_prop*100:.1f}%")
print(f"  Daylight duration:        {day_length*24:.1f} hours\n")

print("=== Weights ===")
print(f"  Payload weight:           {W_payload:.2f} N")
print(f"  Spar weight:              {W_spar_sol:.2f} N")
print(f"  Skin weight:              {sol.value(W_skin):.2f} N")
print(f"  Battery mass:             {battery_mass_sol:.2f} kg")
print(f"  Battery weight:           {sol.value(W_battery):.2f} N")
print(f"  Total weight:             {Total_weight_sol:.2f} N\n")

print("=== Structural ===")
print(f"  Spar thicknesses (main):\n    {t_list_main_sol}")
print(f"  Spar thicknesses (second):\n    {t_list_second_sol}")

print("\n" + "="*35)
