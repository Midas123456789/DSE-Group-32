import aerosandbox as asb
import aerosandbox.numpy as np
from mass_wing import Mass_wing


wing_airfoil = asb.Airfoil("sd7037")
N_cords = 3
g = 9.81
payload_mass = 2000  #kg
Wp = payload_mass * g
altitude = 15000
atm = asb.Atmosphere(altitude=altitude)
opti = asb.Opti()


cords = opti.variable(init_guess=8 * np.ones(N_cords), n_vars=N_cords)
b = opti.variable(init_guess=40, upper_bound=50, lower_bound=0)
V = opti.variable(init_guess=30, upper_bound=200, lower_bound=5)
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

L = 0.5 * atm.density() * V**2 * aero['CL'] * wing.area()
root_r = wing.root_R()
tip_r = wing.tip_R()
spar_mass = wing.spar_mass(L, b)
spar_thickness = wing.root_t(L, b)
Ww = spar_mass * g
W = Wp + Ww
opti.subject_to(
    L > W
)

Pg = wing.area() * 0.2 * 1000
Pr = 0.5 * atm.density() * V**2 * wing.area() * aero['CD'] * V

opti.subject_to(
    Pg > Pr
)

# Objective: Minimize drag power (drag force * velocity)
opti.minimize(
    Pr
)

sol = opti.solve()

# Print results
print("Optimized span b (m):", sol(b))
print("Optimized alpha (deg):", sol(alpha))
print("Optimized chords (m):", sol(cords))
print("Wing area:", sol(wing.area()))
print("CL:", sol(aero['CL']))
print("CD:", sol(aero['CD']))
print('V = ', sol(V))
print('root r = ', sol(root_r))
print('tip r = ', sol(tip_r))
print('spar mass = ', sol(spar_mass))
print('spar thickness = ', sol(spar_thickness))
print('W = ', sol(W))
print('L = ', sol(L))
print('Pr = ', sol(Pr))
print('Pg = ', sol(Pg))


# Visualize
vlm = sol(vlm)
vlm.draw(show_kwargs=dict(jupyter_backend="static"))

