import aerosandbox as asb
import aerosandbox.numpy as np
import aerosandbox.tools.units as u
#from mass_wing import Mass_wing

def wing_mass(M_TO, S, n_ult, A):
    """
    Estimate wing weight using the empirical formula:
    
    W_w = 0.04674 * (W_TO)^0.397 * (S)^0.360 * (n_ult)^0.397 * (A)^1.712
    
    Parameters:
        W_TO (float): Take-off weight
        S (float): Wing area
        n_ult (float): Ultimate load factor
        A (float): Aspect ratio
    
    Returns:
        float: Estimated wing weight
    """
    return (0.04674 * ((M_TO / u.lbm) ** 0.397) * ((S / (u.foot ** 2)) ** 0.360) * ((n_ult) ** 0.397) * ((A) ** 1.712)) * u.lbm


# === Constants === #
wing_airfoil = asb.Airfoil("sd7037")
N_cords = 2
g = 9.81
altitude = 15000
atm = asb.Atmosphere(altitude=altitude)
opti = asb.Opti()

# === Solar & Power Parameters === #
solar_irradiance = 700         # W/m^2 available at altitude
solar_efficiency = 0.20        # 20% cell efficiency
motor_eff = 0.8
prop_eff = 0.8
η_prop = motor_eff * prop_eff  # Total propulsive efficiency

# === Mass Parameters === #
W_payload = 1000  # N

# === Design Variables === #
cords = opti.variable(init_guess=3 * np.ones(N_cords), n_vars=N_cords)
b = opti.variable(init_guess=60, upper_bound=70, lower_bound=0)
V = opti.variable(init_guess=30, upper_bound=200, lower_bound=1)
alpha = opti.variable(init_guess=5, lower_bound=0, upper_bound=30)
W_TO = opti.variable(init_guess=2*W_payload, upper_bound=10000, lower_bound=1)

# === Geometry Setup === #
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
).translate([4, 0, 0])
airplane = asb.Airplane(wings=[wing])

# === Constraints === #
opti.subject_to([
    cords > 0,
    np.diff(cords) <= 0,  # Taper
    wing.area() < 2000
])

# === Aerodynamics === #
op = asb.OperatingPoint(velocity=V, alpha=alpha, atmosphere=atm)
aero = asb.AeroBuildup(airplane=airplane, op_point=op).run()

# === Lift & Drag Forces === #
L = 0.5 * atm.density() * V**2 * aero['CL'] * wing.area()
D = 0.5 * atm.density() * V**2 * aero['CD'] * wing.area()

# === Skin Weight === #
skin_density = 0.25  # kg/m², adjust based on material


# === Wing Weight === #
M_TO = W_TO / 9.81
M_wing = wing_mass(M_TO, wing.area(), 3, wing.aspect_ratio()) 
W_wing = g * M_wing
# === Total Weight === #
opti.subject_to(L == W_TO)
opti.subject_to(W_TO == W_wing + W_payload)

# === Payload Power === #
P_payload = 5000  # Watts, constant power for onboard systems

# === Power Required & Generated === #
P_required = D * V
P_generated = wing.area() * solar_irradiance * solar_efficiency
P_available = P_generated * η_prop - P_payload  # Account for payload power draw

# === Power Constraint === #
opti.subject_to(P_available >= P_required)

# === Objective === #
opti.minimize(P_required)

# === Solve === #
sol = opti.solve()

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
P_generated_sol = sol.value(P_generated)
P_available_sol = sol.value(P_available)
Wing_weight_sol = sol.value(W_wing)
Total_weight_sol = sol.value(W_TO)
AR_sol = sol.value(wing.aspect_ratio())

# === Rebuild Solved Wing === #
y_sections_sol = np.linspace(0, b_sol / 2, N_cords)
wing_sol = asb.Wing(
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

airplane_sol = asb.Airplane(wings=[wing_sol])
airplane_sol.draw(show=True)

# === Print Results === #
print("==== Optimization Results ====")
print(f"Span b:               {b_sol:.2f} m")
print(f"Chord distribution:   {cords_sol}")
print(f"Flight velocity:      {V_sol:.2f} m/s")
print(f"Angle of attack:      {alpha_sol:.2f} deg")
print(f"Wing area:            {area_sol:.2f} m²")
print('AR: ', AR_sol)
print(f"CL:                   {CL_sol:.3f}")
print(f"CD:                   {CD_sol:.3f}")
print(f"Lift:                 {L_sol:.2f} N")
print(f"Drag:                 {D_sol:.2f} N")
print(f"Power required:       {P_required_sol:.2f} W")
print(f"Power generated:      {P_generated_sol:.2f} W")
print(f"Power available:      {P_available_sol:.2f} W")
print('Wing weight = ', Wing_weight_sol)
print('Total weight = ', Total_weight_sol)


