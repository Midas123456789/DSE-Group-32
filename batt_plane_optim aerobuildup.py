import aerosandbox as asb
import aerosandbox.numpy as np
import matplotlib.pyplot as plt
from mass_wing import Mass_wing

# === Constants === #
g = 9.81  # m/s²
altitude = 15000  # m
payload_mass = 500  # kg
Wp = payload_mass * g  # Payload weight in N
ribs = 10
skin_density = 0.220  # kg/m²
Ppl = 5000  # Payload power in W

# Battery properties
battery_days = 300
battery_cost = 20  # $/kg
bat_energy_density = 576000  # J/kg

# Efficiencies
propeller_eff = 0.8
motor_eff = 0.8
prop_eff = propeller_eff * motor_eff

# Atmosphere
atm = asb.Atmosphere(altitude=altitude)

# === Optimization Setup === #
opti = asb.Opti()
N_cords = 3

cords = opti.variable(init_guess=3 * np.ones(N_cords), n_vars=N_cords)
b = opti.variable(init_guess=60, upper_bound=70, lower_bound=0)  # Initial guess < upper bound
V = opti.variable(init_guess=30, upper_bound=200, lower_bound=40)
alpha = opti.variable(init_guess=5, lower_bound=0, upper_bound=30)

# Spanwise locations
y_sections = np.linspace(0, b / 2, N_cords)

# Airfoil and Wing
wing_airfoil = asb.Airfoil("sd7037")
wing = Mass_wing(
    symmetric=True,
    xsecs=[
        asb.WingXSec(
            xyz_le=[-0.25 * cords[i], y_sections[i], 0],
            chord=cords[i],
            airfoil=wing_airfoil
        ) for i in range(N_cords)
    ]
).translate([4, 0, 0])

airplane = asb.Airplane(wings=[wing])

# === Constraints === #
opti.subject_to([
    cords > 0,
    np.diff(cords) <= 0,  # Taper constraint
    wing.area() < 2000
])

# Operating point and aerodynamics
op = asb.OperatingPoint(velocity=V, alpha=alpha, atmosphere=atm)
aero = asb.AeroBuildup(airplane=airplane, op_point=op).run()

# Power and Energy
Pg = wing.area() * 0.2 * 700 * 0.5  # 20% efficient panels, 700 W/m² solar, 50% day
Pr = 0.5 * atm.density() * V**2 * wing.area() * aero['CD'] * V + Ppl
night_energy = Pr * 0.5 * 86400  # 12-hour night
Mbat = night_energy / bat_energy_density
Wbat = Mbat * g

# Structural Masses
Mskin = 2 * wing.area() * skin_density
Wskin = Mskin * g
spar_mass = wing.spar_mass(L := 0.5 * atm.density() * V**2 * aero['CL'] * wing.area(), b)
Ww = spar_mass * g

# Total weight
W = Wp + Ww + Wbat + Wskin

# Final constraints
opti.subject_to([
    Pg * prop_eff > Pr,  # Solar power must exceed required power
    L > W                # Lift must exceed weight
])

# Objective
opti.minimize(W)

# === Solve === #
sol = opti.solve()

# === Post-processing === #
cords_sol = sol.value(cords)
b_sol = sol.value(b)
V_sol = sol.value(V)
alpha_sol = sol.value(alpha)
area_sol = sol.value(wing.area())
CL_sol = sol.value(aero['CL'])
CD_sol = sol.value(aero['CD'])
L_sol = sol.value(L)
Pr_sol = sol.value(Pr)

# Rebuild optimized wing
y_sections_sol = np.linspace(0, b_sol / 2, N_cords)
wing_sol = Mass_wing(
    symmetric=True,
    xsecs=[
        asb.WingXSec(
            xyz_le=[-0.25 * cords_sol[i], y_sections_sol[i], 0],
            chord=cords_sol[i],
            airfoil=wing_airfoil
        ) for i in range(N_cords)
    ]
).translate([4, 0, 0])
airplane_sol = asb.Airplane(wings=[wing_sol])
airplane_sol.draw(show=True)

# === Print Results === #
print("==== Optimization Results ====")
print(f"Span b:               {b_sol:.2f} m")
print(f"Chord distribution:   {cords_sol}")
print(f"Flight velocity:      {V_sol:.2f} m/s")
print(f"Angle of attack:      {alpha_sol:.2f} deg")
print(f"Wing area:            {area_sol:.2f} m²")
print(f"CL:                   {CL_sol:.3f}")
print(f"CD:                   {CD_sol:.3f}")
print(f"Lift:                 {L_sol:.2f} N")
print(f"Power required:       {Pr_sol:.2f} W")

# Optional breakdown
print("\n---- Mass Breakdown ----")
print(f"Payload weight:       {Wp:.2f} N")
print(f"Wing structure:       {sol.value(Ww):.2f} N")
print(f"Battery weight:       {sol.value(Wbat):.2f} N")
print(f"Skin weight:          {sol.value(Wskin):.2f} N")
print(f"Total weight:         {sol.value(W):.2f} N")

# === Plot Chord Distribution === #
plt.figure()
plt.plot(y_sections_sol, cords_sol, marker="o")
plt.xlabel("Spanwise position y [m]")
plt.ylabel("Chord length [m]")
plt.title("Optimized Chord Distribution")
plt.grid(True)
plt.tight_layout()
plt.show()
