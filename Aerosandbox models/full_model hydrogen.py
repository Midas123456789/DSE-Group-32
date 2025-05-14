import aerosandbox as asb
import aerosandbox.numpy as np
from mass_wing import Mass_wing

# === Constants === #
wing_airfoil = asb.Airfoil("sd7037")
N_cords = 5
g = 9.81
altitude = 18000
atm = asb.Atmosphere(altitude=altitude)
opti = asb.Opti()


Tail_drag_factor = 1.1

second_in_day = 86400
mission_days = 10
mission_seconds = mission_days * second_in_day

W_payload = 1000  # N
P_payload = 10000  # Watts, constant power for onboard systems

density_LH = 70.85 #kg/m3
energy_density_LH = 142 * 10 ** 6 #J/kg
#fuel_tank_wall_thickness = 0.0612  # from Brewer, Hydrogen Aircraft Technology pg. 203
power_density_PEMFCs = 40000 #W/kg
PEMFCs_eff = 0.55
motor_eff = 0.97 #Fundamentals of Aircraft and Airship Design
propeller_eff = 0.85 #Fundamentals of Aircraft and Airship Design
conversion_eff = PEMFCs_eff
propulsive_eff = conversion_eff * motor_eff * propeller_eff
fuel_tank_fuel_mass_fraction = 0.34  # from Brewer, Hydrogen Aircraft Technology pg. 30

# === Design Variables === #
cords = opti.variable(init_guess=3 * np.ones(N_cords), n_vars=N_cords)
b = opti.variable(init_guess=60, upper_bound=200, lower_bound=0)
V = opti.variable(init_guess=30, upper_bound=200, lower_bound=1)
alpha = opti.variable(init_guess=5, lower_bound=0, upper_bound=30)

# === Geometry Setup === #
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
).translate([0, 0, 0])

tail_scalar = 1
def scale(vec, scalar):
    return [x * scalar for x in vec]

tail = [
    asb.Wing(
        name="H-stab",
        symmetric=True,
        xsecs=[
            asb.WingXSec(
                xyz_le=scale([0, 0, 0], tail_scalar),
                chord=tail_scalar * 2,
                airfoil=asb.Airfoil("ht08")
            ),
            asb.WingXSec(
                xyz_le=scale([0, 3, 0], tail_scalar),
                chord=tail_scalar * 1,
                airfoil=asb.Airfoil("ht08")
            ),
        ]
    ).translate([15, 0, 0]),

    asb.Wing(
        name="V-stab",
        xsecs=[
            asb.WingXSec(
                xyz_le=scale([0, 0, 0], tail_scalar),
                chord=tail_scalar * 2,
                airfoil=asb.Airfoil("ht08")
            ),
            asb.WingXSec(
                xyz_le=scale([0.5, 0, 3], tail_scalar),
                chord=tail_scalar * 1,
                airfoil=asb.Airfoil("ht08")
            )
        ]
    ).translate([16, 0, 0])
]
# === Constraints === #
opti.subject_to([
    cords > 0,
    np.diff(cords) <= 0,  # Taper
    wing.area() < 2000
])
if False:
    wings = [wing] + tail
else:
    wings = [wing]
airplane = asb.Airplane(wings=wings)


# === Aerodynamics === #
op = asb.OperatingPoint(velocity=V, alpha=alpha, atmosphere=atm)
aero = asb.AeroBuildup(airplane=airplane, op_point=op).run()

# === Lift & Drag Forces === #
L = 0.5 * atm.density() * V**2 * aero['CL'] * wing.area()
D = 0.5 * atm.density() * V**2 * aero['CD'] * wing.area() * Tail_drag_factor

# === Payload Power === #

# === Power Constraint and Battery Sizing === #
P_required = V * D
E_LH = 1 / conversion_eff * (P_required / propulsive_eff + P_payload) * mission_seconds
M_LH = E_LH / energy_density_LH
W_LH = M_LH * g
M_tank = fuel_tank_fuel_mass_fraction * M_LH
W_tank = M_tank * g
LH_volume = M_LH / density_LH

# === Skin Weight === #
skin_density = 0.25  # kg/m², adjust based on material
W_skin = 2 * wing.area() * skin_density * g

# === Wing Weight === #
spar_mass, t_list = wing.spar_mass(L, b)
W_spar = spar_mass * g   
W_wing = W_spar + W_skin
M_fuel_cell = E_LH / mission_seconds / power_density_PEMFCs
W_fuel_cell = M_fuel_cell * g
# === Miss Weight === #
W_misc = W_wing

# === Total Weight === #
W_total = W_payload + W_wing + W_LH + W_misc + W_tank + W_fuel_cell
opti.subject_to(L >= W_total)

# === Objective === #
opti.minimize(P_required)
# === Solve === #
sol = opti.solve(verbose=False)

# === Extract Solved Values === #
cords_sol = sol.value(cords)
b_sol = sol.value(b)
V_sol = sol.value(V)
alpha_sol = sol.value(alpha)
area_sol = sol.value(wing.area())
CL_sol = sol.value(aero['CL'])
CD_sol = sol.value(aero['CD'])
L_sol = sol.value(L)
D_sol = sol.value(D)
P_required_sol = sol.value(P_required)
#P_generated_sol = sol.value(P_generated)
#P_available_sol = sol.value(P_available)
W_spar_sol = sol.value(W_spar)
t_list_sol = sol(t_list)
Total_weight_sol = sol(W_total)
W_misc = sol(W_misc)

# === Rebuild Solved Wing === #
y_sections_sol = np.linspace(0, b_sol / 2, N_cords)
wing_sol = Mass_wing(
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

airplane_sol = asb.Airplane(wings=[wing_sol] + tail)
# === Output Results === #
print("\n" + "="*35)
print("        OPTIMIZATION RESULTS")
print("="*35 + "\n")

print("=== Geometry ===")
print(f"  Span (b):                 {b_sol:.2f} m")
print(f"  Chord distribution (m):   {cords_sol}")
print(f"  Wing area:                {area_sol:.2f} m²")
print(f"  Volume LH:                {sol.value(LH_volume):.2f} m3\n")

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
print(f"  Drag (L/D):                 {L_sol/D_sol:.2f} \n")

print("=== Energy & Power ===")
print(f"  Payload power:            {P_payload:.2f} W")
print(f"  Power required (P):       {P_required_sol:.2f} W")


print("=== Weights ===")
print(f"  Payload weight:           {W_payload:.2f} N")
print(f"  Spar weight:              {W_spar_sol:.2f} N")
print(f"  Skin weight:              {sol.value(W_skin):.2f} N")
print(f"  LH weight:           {sol.value(W_LH):.2f} N")
print(f"  Miscellaneous weight:     {W_misc:.2f} N")  # New line for miscellaneous weight
print(f"  Total weight:             {Total_weight_sol:.2f} N\n")

print("=== Mass ===")
print(f"  Payload mass:           {W_payload / 9.81:.2f} kg")
print(f"  Spar mass:              {W_spar_sol / 9.81:.2f} kg")
print(f"  Skin mass:              {sol.value(W_skin) / 9.81:.2f} kg")
print(f"  Miscellaneous mass:     {W_misc / 9.81:.2f} kg")  # New line for miscellaneous weight
print(f"  LH mass:           {sol.value(M_LH):.2f} kg")
print(f"  Fuel cell mass:           {sol.value(M_fuel_cell):.2f} kg")
print(f"  Tank mass:           {sol.value(M_tank):.2f} kg")
print(f"  Total mass:             {Total_weight_sol / 9.81:.2f} kg\n")

print('=== Operations ===')
print(f"  Yearly LH consumption:           {sol.value(M_LH)*52:.2f} kg")
eu_per_kg_LH = 10
print(f"  Yearly Fuel Cost:           {sol.value(M_LH)*52 * eu_per_kg_LH:.2f} Euro\n")

print("=== Structural ===")
print(f"  Spar thicknesses (t_list):\n    {t_list_sol}")

print("\n" + "="*35)
airplane_sol.draw(show=True)