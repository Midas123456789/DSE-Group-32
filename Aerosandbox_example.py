import aerosandbox as asb
import aerosandbox.numpy as np

wing_airfoil = asb.Airfoil("sd7037")
N_cords = 4
g = 9.81
payload_mass = 3962  #kg
W = payload_mass * g
V = 30  # m/s (constant for now)
solar_intensity = 1380
solar_eff_fact = 0.2
altitude = 20000
atm = asb.Atmosphere(altitude=altitude)

opti = asb.Opti()

cords = opti.variable(init_guess=8 * np.ones(N_cords), n_vars=N_cords)
b = opti.variable(init_guess=40, upper_bound=50, lower_bound=30)

# Define y-locations based on variable b
y_sections = np.linspace(0, b / 2, N_cords)

wing = asb.Wing(
    symmetric=True,
    xsecs=[
        asb.WingXSec(
            xyz_le=[-0.25 * cords[i], y_sections[i], 0],
            chord=cords[i],
            airfoil=wing_airfoil
        )
        for i in range(N_cords)
    ]
)

airplane = asb.Airplane(wings=[wing])

opti.subject_to([
    cords > 0,
    wing.area() < 200,
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

k = 0.1 # Structural weight coefficient [N/m³]
wing_weight = k * b**2

# Lift = Total Weight
opti.subject_to(
    0.5 * atm.density() * V**2 * aero['CL'] * wing.area() == W
)


opti.subject_to(
    0.5 * atm.density() * V**2 * wing.area() * aero['CD'] * V < solar_eff_fact * solar_intensity * wing.area()
)

opti.minimize(
    b
)

sol = opti.solve()

# Print results
print("Optimized span b (m):", sol(b))
print("Optimized alpha (deg):", sol(alpha))
print("Optimized chords (m):", sol(cords))
print("Wing area:", sol(wing.area()))
print("CL:", sol(aero['CL']))
print("CD:", sol(aero['CD']))
print("Wing structural weight (N):", k * sol(b) * sol(wing.area()))
print("Power available: ", solar_eff_fact * solar_intensity * sol(wing.area()))
print("Power req: ", 0.5 * atm.density() * V**2 * sol(wing.area()) * sol(aero['CD']) * V)

# Visualize
vlm = sol(vlm)
vlm.draw(show_kwargs=dict(jupyter_backend="static"))
