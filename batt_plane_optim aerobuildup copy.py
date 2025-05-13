import aerosandbox as asb
import aerosandbox.numpy as np  # Symbolic np for optimization

opti = asb.Opti()

### --- Variables ---
W_dg = opti.variable(init_guess=1000)  # Design gross weight [N]
S_w = opti.variable(init_guess=20)     # Wing area [m^2]
A = opti.variable(init_guess=8)        # Aspect ratio
t_c_root = opti.variable(init_guess=0.12)
taper_ratio = opti.variable(init_guess=0.5)
sweep_deg = 0
sweep_rad = 0
S_csw = opti.variable(init_guess=1.0)

# Assume load factor is 1 for straight and level flight
N_z = 1.0
rho = 1.225  # Air density [kg/m^3]
V = 30       # Flight speed [m/s]

### --- Custom Wing Weight Formula ---
W_wing = (
    0.0051 *
    (W_dg + N_z) ** 0.557 *
    S_w ** 0.649 *
    A ** 0.5 *
    t_c_root ** -0.4 *
    (1 + taper_ratio) ** 0.1 *
    (np.cos(sweep_rad)) ** -1.0 *
    S_csw ** 0.1
)

### --- Geometry Setup ---
wing = asb.Wing(
    name="Main Wing",
    xsecs=[
        asb.WingXSec(
            xyz_le=[0, 0, 0],
            chord=2,
            twist=0,
            airfoil=asb.Airfoil("naca2412"),
        ),
        asb.WingXSec(
            xyz_le=[4, 5, 0],  # span-wise station
            chord=1,
            twist=0,
            airfoil=asb.Airfoil("naca2412"),
        ),
    ]
)

### --- Aerodynamics ---
aero = asb.aeroBuildup(
    airplane=asb.Airplane(
        wings=[wing],
        mass_props=asb.MassProperties(
            weight=W_dg,
        ),
    ),
    op_point=asb.OperatingPoint(
        velocity=V,
        density=rho,
        alpha=opti.variable(init_guess=5),
    )
)

### --- Lift equals Weight constraint ---
opti.subject_to(aero["L"] == W_dg)

### --- Objective ---
opti.minimize(W_wing)  # Example: minimize wing weight

### --- Solve ---
sol = opti.solve()
print("Optimal Wing Weight:", sol.value(W_wing))
print("Optimal Design Gross Weight:", sol.value(W_dg))
print("Lift generated:", sol.value(aero["L"]))
