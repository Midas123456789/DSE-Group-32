import aerosandbox as asb
import aerosandbox.numpy as np
from aerosandbox import cas
from numpy import pi, sqrt
from aerosandbox import Atmosphere
# Constant
altitude=20000

D = 300       # Thrust [N]
v_o=40        # m/s




atmo = Atmosphere(altitude)

rho = atmo.density() # Air density [kg/m^3]

# AeroSandbox Optimization
opti = asb.Opti()

# Optimization variables

r = opti.variable(init_guess=1, lower_bound=0.1, upper_bound=2)    # m

# Derived values
Adisk = pi * r**2

# Power required (momentum theory)
T = D
power = 0.5 * rho * T * v_o * (cas.sqrt(T / (0.5 * rho * Adisk * v_o**2) + 1) + 1)

# Objective: Minimize power
opti.minimize(power)

# Solve
sol = opti.solve()

# Get optimized values

r_opt = sol.value(r)
P_opt = sol.value(power)

# Output results
print(f"Optimized flight speed (v_o): {v_o:.2f} m/s")
print(f"Optimized propeller radius (r): {r_opt:.2f} m")
print(f"Minimum power required: {P_opt:.2f} W")
print("density: ", rho)