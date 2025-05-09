import aerosandbox as asb
import aerosandbox.numpy as np
import matplotlib.pyplot as plt
from aerosandbox import Atmosphere

# Inputs
#Design variables
W = 7000  # Weight [N]
D = 0.05*W    # Drag [N]
V = 40     # Speed [m/s]

#Mission variables
h_target = 18000  # Target altitude [m]
night_time = 12  # [hr]
day_time = 24 - night_time  # [hr]
mission_days = 30
mission_hours = mission_days * 24  # [hr]
energy_margin = 1.05

#Energy parameters
solar_irradiance = 800  # [W/m^2]
panel_efficiency = 0.2
energy_density_bat=450  # [Wh/kg]

#Legal limits
Maximum_possible_descend_altitude=14000

#CODE


opti = asb.Opti()

T = opti.variable(init_guess=1.4 * D, lower_bound=D, upper_bound=0.5 * W)
A = opti.variable(init_guess=50, lower_bound=0, upper_bound=600)
b_capacity = opti.variable(init_guess=80000, lower_bound=0, upper_bound=500000)
h_difference = opti.variable(init_guess=1000, lower_bound=0, upper_bound=(h_target - Maximum_possible_descend_altitude))
y = opti.variable(init_guess=0.1, lower_bound=0.01, upper_bound=0.2)

# Climb
RC = (T * V - D * V) / W
t_climb = h_target / RC  # [s]
E_climb = T * V * t_climb  # [J]

# Descent
E_descent_gain = h_difference * W  # [J]

# Cruise
E_cruise_day = D * V * (day_time * 3600 - t_climb)  # [J]

# Day total
E_day_total = (E_climb + E_cruise_day ) * energy_margin
E_day_total_Wh = E_day_total / 3600

# Night energy required
t_descent = h_difference / (V * np.sin(y))
t_cruise_night = night_time * 3600 - t_descent
E_night_cruise = D * V * t_cruise_night  # J
E_night_net = E_night_cruise - h_difference * W  # Net required energy
E_night_Wh = (E_night_net * energy_margin) / 3600

# Solar energy input
solar_input_Wh = solar_irradiance * A * panel_efficiency * day_time

# Constraints
opti.subject_to(b_capacity >= E_night_Wh)
opti.subject_to(solar_input_Wh >= E_day_total_Wh + E_night_Wh)

# Objective
opti.minimize(A)


sol = opti.solve()

# Extract
b_capacity_opt = sol.value(b_capacity)
solar_area_opt = sol.value(A)
thrust_opt = sol.value(T)
h_difference_opt = sol.value(h_difference)
y_opt = sol.value(y)
RC_opt = (thrust_opt * V - D * V) / W
t_climb_opt = h_target / RC_opt
t_descent_opt = h_difference_opt / (V * np.sin(y_opt))

# Altitude profile over 72 hrs
altitude_profile = []
time_profile = []
alt = 0
time_hr = 0

for day in range(mission_hours // 24):
    ### DAY ###
    # Climb
    t_climb_hr = t_climb_opt / 3600
    climb_rate = h_target / t_climb_hr
    t_vals = np.linspace(0, t_climb_hr, 20)
    h_vals = np.linspace(alt, h_target, 20)
    time_profile.extend(time_hr + t_vals)
    altitude_profile.extend(h_vals)
    time_hr += t_climb_hr
    alt = h_target

    # Cruise
    t_cruise = day_time - t_climb_hr
    t_vals = np.linspace(0, t_cruise, 10)
    h_vals = np.full_like(t_vals, alt)
    time_profile.extend(time_hr + t_vals)
    altitude_profile.extend(h_vals)
    time_hr += t_cruise

    ### NIGHT ###
    # Descent
    t_descent_hr = t_descent_opt / 3600
    descent_target = h_target - h_difference_opt
    t_vals = np.linspace(0, t_descent_hr, 20)
    h_vals = np.linspace(alt, descent_target, 20)
    time_profile.extend(time_hr + t_vals)
    altitude_profile.extend(h_vals)
    time_hr += t_descent_hr
    alt = descent_target

    # Cruise at lower altitude
    t_cruise_night = night_time - t_descent_hr
    t_vals = np.linspace(0, t_cruise_night, 10)
    h_vals = np.full_like(t_vals, alt)
    time_profile.extend(time_hr + t_vals)
    altitude_profile.extend(h_vals)
    time_hr += t_cruise_night

final_descent_time_hr = alt / (V * np.sin(y_opt)) / 3600  # time from current altitude to 0
t_vals = np.linspace(0, final_descent_time_hr, 30)
h_vals = np.linspace(alt, 0, 30)
time_profile.extend(time_hr + t_vals)
altitude_profile.extend(h_vals)
time_hr += final_descent_time_hr

# Plot
plt.figure(figsize=(12, 6))
plt.plot(time_profile, altitude_profile)
plt.title(f"{mission_hours}-Hour Mission Profile: Altitude vs Time")
plt.xlabel("Time [hr]")
plt.ylabel("Altitude [m]")
plt.grid(True)
plt.tight_layout()
plt.show()

# Summary
print("Battery capacity (Wh):", round(b_capacity_opt, 1))
print("Battery weight (kg):", round(b_capacity_opt/energy_density_bat, 1))
print("Battery weight percentage of MTOW %:", round(((b_capacity_opt/energy_density_bat)/(W/9.81))*100, 3))
print("Solar panel area (m^2):", round(solar_area_opt, 1))
print("Thrust (N):", round(thrust_opt, 1))
print("Altitude difference during night descent (m):", round(h_difference_opt, 1))
print("Decent angle (deg):", round(y_opt*180/np.pi, 4))
print("Rate of climb (m/s):", round(RC_opt, 2))
