import aerosandbox as asb
import aerosandbox.numpy as np
from mass_wing import Mass_wing



wing_airfoil = asb.Airfoil("sd7037")
N_cords = 5
g = 9.81
payload_mass = 500  #kg
Wp = payload_mass * g
altitude = 15000
ribs = 10
skin_density = 0.220 #kg/m2
atm = asb.Atmosphere(altitude=altitude)
opti = asb.Opti()


cords = opti.variable(init_guess=8 * np.ones(N_cords), n_vars=N_cords)
b = opti.variable(init_guess=40, upper_bound=50, lower_bound=0)
V = opti.variable(init_guess=30, upper_bound=200, lower_bound=40)
# Define y-locations based on variable b
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
    wing.area() < 2000,
    np.diff(cords) <= 0,  # Enforce taper
])

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

Pg = wing.area() * 0.2 * 500 * 0.5
Pr = 0.5 * atm.density() * V**2 * wing.area() * aero['CD'] * V
bat_energy_density = 576000
night_energy = Pr * 0.5 * 86400
Mbat = night_energy / bat_energy_density
Wbat = Mbat * g

Mskin = 2 * wing.area() * skin_density
Wskin = Mskin * g
opti.subject_to(
    Pg > Pr
)

# Objective: Minimize drag power (drag force * velocity)
opti.minimize(
    Pr
)


L = 0.5 * atm.density() * V**2 * aero['CL'] * wing.area()
root_r = wing.root_R()
tip_r = wing.tip_R()
spar_mass = wing.spar_mass(L, b)
spar_thickness = wing.root_t(L, b)
Ww = spar_mass * g



W = Wp + Ww + Wbat + Wskin
opti.subject_to(
    L > W
)


sol = opti.solve()

# Print results
# Geometry
print("=== GEOMETRY ===")
print("Optimized span b (m):", sol(b))
print("Optimized alpha (deg):", sol(alpha))
print("Optimized chords (m):", sol(cords))
print("Wing area:", sol(wing.area()))
print("Root r =", sol(root_r))
print("Tip r =", sol(tip_r))
print()

# Aerodynamics
print("=== AERODYNAMICS ===")
print("CL:", sol(aero['CL']))
print("CD:", sol(aero['CD']))
print()

# Performance
print("=== PERFORMANCE ===")
print("V =", sol(V))
print("Lift L =", sol(L))
print("Required Power Pr =", sol(Pr))
print("Generated Power Pg =", sol(Pg))
print()

# Structure
print("=== STRUCTURE ===")
print("Spar thickness =", sol(spar_thickness))

print("=== MASS ===")
print("Mass =", sol(W) / g)
print("Spar mass =", sol(spar_mass))
print("Mass bat =", sol(Mbat))
print("Mass skin =", sol(Mskin))
print("Wing loading Kg/m2 =", sol(W) / g / sol(wing.area()))


# Visualize
vlm = sol(vlm)
vlm.draw(show_kwargs=dict(jupyter_backend="static"))

