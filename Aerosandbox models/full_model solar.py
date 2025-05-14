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

# === Solar & Power Parameters === #
solar_irradiance = 120         # W/m^2 available at altitude
solar_efficiency = 0.30        # 20% cell efficiency
motor_eff = 0.97 #Fundamentals of Aircraft and Airship Design
prop_eff = 0.85 #Fundamentals of Aircraft and Airship Design
η_prop = motor_eff * prop_eff  # Total propulsive efficiency
Tail_drag_factor = 1.1


W_payload = 500  # N
P_payload = 2e3 # Watts, constant power for onboard systems

Energy_density_bat = 1.8e6  # J/kg # based on Elysian

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
    ).translate([10, 0, 0]),

    asb.Wing(
        name="V-stab",
        xsecs=[
            asb.WingXSec(
                xyz_le=scale([0, 0, 0], tail_scalar),
                chord=tail_scalar * 2,
                airfoil=asb.Airfoil("ht08")
            ),
            asb.WingXSec(
                xyz_le=scale([0.5, 0, 2], tail_scalar),
                chord=tail_scalar * 1,
                airfoil=asb.Airfoil("ht08")
            )
        ]
    ).translate([10, 0, 0])
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
day_frac = 0.5
second_in_day = 86400
P_required = V * D
E_generated = wing.area() * solar_irradiance * solar_efficiency * day_frac * second_in_day
E_payload = P_payload * second_in_day
E_prop = P_required * second_in_day / η_prop
E_bat = (E_prop + E_payload) * (1-day_frac)
M_bat = E_bat / Energy_density_bat
W_battery = M_bat * g

opti.subject_to(
    E_generated > E_payload + E_prop
)

# === Skin Weight === #
skin_density = 0.25  # kg/m², adjust based on material
W_skin = 2 * wing.area() * skin_density * g

# === Wing Weight === #
spar_mass, t_list = wing.spar_mass(L, b)
W_spar = spar_mass * g   
W_wing = W_spar + W_skin

# === Miss Weight === #
W_misc = W_wing

# === Total Weight === #
W_total = W_payload + W_wing + W_battery + W_misc
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
battery_mass_sol = sol(M_bat)
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
print(f"  Wing area:                {area_sol:.2f} m²\n")

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

print(f"  Solar efficiency:         {solar_efficiency*100:.1f}%")
print(f"  Propulsion efficiency:    {η_prop*100:.1f}%")
print(f"  Daylight duration:        {day_frac*24:.1f} hours\n")

print("=== Weights ===")
print(f"  Payload weight:           {W_payload:.2f} N")
print(f"  Spar weight:              {W_spar_sol:.2f} N")
print(f"  Skin weight:              {sol.value(W_skin):.2f} N")
print(f"  Battery weight:           {sol.value(W_battery):.2f} N")
print(f"  Miscellaneous weight:     {W_misc:.2f} N")  # New line for miscellaneous weight
print(f"  Total weight:             {Total_weight_sol:.2f} N\n")

print("=== Mass ===")
print(f"  Payload mass:           {W_payload / 9.81:.2f} kg")
print(f"  Spar mass:              {W_spar_sol / 9.81:.2f} kg")
print(f"  Skin mass:              {sol.value(W_skin) / 9.81:.2f} kg")
print(f"  Battery mass:             {battery_mass_sol:.2f} kg")
print(f"  Miscellaneous mass:     {W_misc / 9.81:.2f} kg")  # New line for miscellaneous weight
print(f"  Total mass:             {Total_weight_sol / 9.81:.2f} kg\n")

print('=== Operations ===')
battery_lifetime = 300 #days
battery_price = 3.11*10**-5 #eu/J
print(f"  Yearly battery consumption:           {sol.value(M_bat)*365/300:.2f} kg")
print(f"  Yearly battery price:           {sol.value(E_bat)*365/300*battery_price:.2f} Euro\n")

print("=== Structural ===")
print(f"  Spar thicknesses (t_list):\n    {t_list_sol}")

print("\n" + "="*35)
#airplane_sol.draw(show=True)