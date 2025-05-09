import aerosandbox as asb
import aerosandbox.numpy as np
from mass_wing import Mass_wing

# Constants
wing_airfoil = asb.Airfoil("sd7037")
N_cords = 2
g = 9.81
payload_mass = 500  # kg
Wp = payload_mass * g
altitude = 1000
ribs = 10
skin_density = 0.220  # kg/m2
atm = asb.Atmosphere(altitude=altitude)
opti = asb.Opti()

# Variables
cords_1 = opti.variable(init_guess=8 * np.ones(N_cords), n_vars=N_cords)
cords_2 = opti.variable(init_guess=8 * np.ones(N_cords), n_vars=N_cords)
b1 = 30
#opti.variable(init_guess=20, upper_bound=25, lower_bound=0)
b2 = 30
#opti.variable(init_guess=20, upper_bound=25, lower_bound=0)
V = 100
#opti.variable(init_guess=30, upper_bound=200, lower_bound=5)

# y-locations
y_sections_1 = np.linspace(0, b1 / 2, N_cords)
y_sections_2 = np.linspace(0, b2 / 2, N_cords)

# Wings
wing1 = Mass_wing(
    symmetric=True,
    xsecs=[
        asb.WingXSec(
            xyz_le=[-0.25 * cords_1[i], y_sections_1[i], 0],
            chord=cords_1[i],
            airfoil=wing_airfoil
        )
        for i in range(N_cords)
    ]
)

wing2 = Mass_wing(
    symmetric=True,
    xsecs=[
        asb.WingXSec(
            xyz_le=[-0.25 * cords_2[i], y_sections_2[i], 0],
            chord=cords_2[i],
            airfoil=wing_airfoil
        )
        for i in range(N_cords)
    ]
).translate([4, 0, 0])

airplane = asb.Airplane(wings=[wing1, wing2])

# Constraints
opti.subject_to([
    cords_1 > 0,
    cords_2 > 0,
    np.diff(cords_1) <= 0,
    np.diff(cords_2) <= 0,
    wing1.area() + wing2.area() < 2000,
])

# Aero conditions
alpha = opti.variable(init_guess=5, lower_bound=0, upper_bound=30)
op_point = asb.OperatingPoint(
    velocity=V,
    alpha=alpha,
    atmosphere=atm
)

vlm = asb.VortexLatticeMethod(
    airplane=airplane,
    op_point=op_point,
    spanwise_resolution=10,
    chordwise_resolution=10
)

aero = vlm.run()

# Power
Pg = (wing1.area() + wing2.area()) * 0.2 * 500 * 0.5
Pr = 0.5 * atm.density() * V**2 * (wing1.area() + wing2.area()) * aero['CD'] * V
bat_energy_density = 576000
night_energy = Pr * 0.5 * 86400
Mbat = night_energy / bat_energy_density
Wbat = Mbat * g

# Mass
Mskin = 2 * (wing1.area() + wing2.area()) * skin_density
Wskin = Mskin * g

opti.subject_to(Pg > Pr)
opti.minimize(Pr)

# Lift
L = 0.5 * atm.density() * V**2 * aero['CL'] * (wing1.area() + wing2.area())

# Structural mass
spar_mass1 = wing1.spar_mass(L / 2, b1)
spar_mass2 = wing2.spar_mass(L / 2, b2)
Ww = (spar_mass1 + spar_mass2) * g

W = Wp + Ww + Wbat + Wskin
opti.subject_to(L > W)

sol = opti.solve()

# Print results
print("=== GEOMETRY ===")
print("Wing 1 span (m):", sol(b1))
print("Wing 2 span (m):", sol(b2))
print("Alpha (deg):", sol(alpha))
print("Wing 1 chords:", sol(cords_1))
print("Wing 2 chords:", sol(cords_2))
print("Total area:", sol(wing1.area() + wing2.area()))
print()

print("=== AERODYNAMICS ===")
print("CL:", sol(aero['CL']))
print("CD:", sol(aero['CD']))
print()

print("=== PERFORMANCE ===")
print("Velocity V =", sol(V))
print("Lift L =", sol(L))
print("Required Power Pr =", sol(Pr))
print("Generated Power Pg =", sol(Pg))
print()

print("=== STRUCTURE ===")
print("Wing 1 spar mass =", sol(spar_mass1))
print("Wing 2 spar mass =", sol(spar_mass2))
print()

print("=== MASS ===")
print("Total mass =", sol(W) / g)
print("Battery mass =", sol(Mbat))
print("Skin mass =", sol(Mskin))

# Visualize
vlm = sol(vlm)
vlm.draw(show_kwargs=dict(jupyter_backend="static"))
