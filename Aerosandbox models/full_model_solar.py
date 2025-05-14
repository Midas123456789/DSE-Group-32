import aerosandbox as asb
import aerosandbox.numpy as np
from mass_wing import Mass_wing

def get_SOL_BAT_conventional(N_cords = 5, wing_airfoil = asb.Airfoil("sd7037"), altitude = 18000, solar_irradiance = 700, W_payload = 1000, P_payload = 10e3):
    g = 9.81
    atm = asb.Atmosphere(altitude=altitude)
    opti = asb.Opti()

    # === Solar & Power Parameters === #
    solar_efficiency = 0.30        # 20% cell efficiency
    motor_eff = 0.97 #Fundamentals of Aircraft and Airship Design
    prop_eff = 0.85 #Fundamentals of Aircraft and Airship Design
    η_prop = motor_eff * prop_eff  # Total propulsive efficiency
    Tail_drag_factor = 1.1

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

    # === Constraints === #
    opti.subject_to([
        cords > 0,
        np.diff(cords) <= 0,  # Taper
        wing.area() < 2000
    ])

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
    
    wing_area = wing.area()
    CL = aero['CL']
    CD = aero['CD'] * Tail_drag_factor
    return {
    "sol": sol,
    "cords": cords,
    "b": b,
    "V": V,
    "alpha": alpha,
    "wing_area": wing_area,
    "CL": CL,
    "CD": CD,
    "L": L,
    "D": D,
    "P_required": P_required,
    "W_spar": W_spar,
    "t_list": t_list,
    "W_total": W_total,
    "W_misc": W_misc,
    "W_skin": W_skin,
    "Airfoil": wing_airfoil,
    "Altitude": altitude,
    "Atm": atm,
    "Solar_panel_area": wing_area,
    "Battery_mass": M_bat,
    "Battery_weight": W_battery
    } 
    
if __name__ == "__main__":
    results = get_SOL_BAT_conventional()
    sol = results["sol"]
    cords_sol = sol.value(results["cords"])
    b_sol = sol.value(results["b"])
    V_sol = sol.value(results["V"])
    alpha_sol = sol.value(results["alpha"])
    area_sol = sol.value(results["wing_area"])
    CL_sol = sol.value(results["CL"])
    CD_sol = sol.value(results["CD"])
    L_sol = sol.value(results["L"])
    D_sol = sol.value(results["D"])
    P_required_sol = sol.value(results["P_required"])
    W_spar_sol = sol.value(results["W_spar"])
    t_list_sol = sol(results["t_list"])
    W_total_sol = sol(results["W_total"])
    W_misc_sol = sol(results["W_misc"])
    W_skin_sol = sol.value(results["W_skin"])
    altitude = sol.value(results["Altitude"])
    atm = sol.value(results["Atm"])
    wing_airfoil = sol.value(results["Airfoil"])
    Solar_area = sol.value(results["Solar_panel_area"])
    Battery_mass = sol.value(results["Battery_mass"])
    Battery_weight = sol.value(results["Battery_weight"])

    # Add manually-defined parameters
    W_payload = 50  # example: payload weight in N
    P_payload = 30  # example: payload power in W
    Total_weight_sol = W_total_sol

    # === Rebuild Solved Wing ===
    y_sections_sol = np.linspace(0, b_sol / 2, len(cords_sol))
    wing_sol = Mass_wing(
        symmetric=True,
        xsecs=[
            asb.WingXSec(
                xyz_le=[-0.25 * cords_sol[i], y_sections_sol[i], 0],
                chord=cords_sol[i],
                airfoil=wing_airfoil
            )
            for i in range(len(cords_sol))
        ]
    ).translate([4, 0, 0])

    airplane_sol = asb.Airplane(wings=[wing_sol])

    # === Output Results ===
    print("\n" + "="*35)
    print("        OPTIMIZATION RESULTS")
    print("="*35 + "\n")

    print("=== Geometry ===")
    print(f"  Span (b):                 {b_sol:.2f} m")
    print(f"  Chord distribution (m):   {cords_sol}")
    print(f"  Wing area:                {area_sol:.2f} m²")
    print(f"  Solar area:                {Solar_area:.2f} m²")

    print("=== Flight Conditions ===")
    print(f"  Flight velocity (V):      {V_sol:.2f} m/s")
    print(f"  Angle of attack (α):      {alpha_sol:.2f}°")
    print(f"  Altitude:                 {altitude} m")
    print(f"  Air density:              {atm.density():.4f} kg/m³\n")

    print("=== Aerodynamics ===")
    print(f"  Lift coefficient (CL):    {CL_sol:.3f}")
    print(f"  Drag coefficient (CD):    {CD_sol:.3f}")
    print(f"  Lift (L):                 {L_sol:.2f} N")
    print(f"  Drag (D):                 {D_sol:.2f} N")
    print(f"  L/D ratio:                {L_sol / D_sol:.2f} \n")

    print("=== Energy & Power ===")
    print(f"  Payload power:            {P_payload:.2f} W")
    print(f"  Power required (P):       {P_required_sol:.2f} W\n")

    print("=== Weights ===")
    print(f"  Payload weight:           {W_payload:.2f} N")
    print(f"  Spar weight:              {W_spar_sol:.2f} N")
    print(f"  Skin weight:              {W_skin_sol:.2f} N")
    print(f"  Miscellaneous weight:     {W_misc_sol:.2f} N")
    print(f"  Battery weight:            {Battery_weight:.2f} N\n")
    print(f"  Total weight:             {Total_weight_sol:.2f} N\n")


    print("=== Mass ===")
    print(f"  Payload mass:             {W_payload / 9.81:.2f} kg")
    print(f"  Spar mass:                {W_spar_sol / 9.81:.2f} kg")
    print(f"  Skin mass:                {W_skin_sol / 9.81:.2f} kg")
    print(f"  Miscellaneous mass:       {W_misc_sol / 9.81:.2f} kg")
    print(f"  Battery mass:       {Battery_mass:.2f} kg")
    print(f"  Total mass:               {Total_weight_sol / 9.81:.2f} kg\n")


    print("=== Structural ===")
    print(f"  Spar thicknesses (t_list):\n    {t_list_sol}")

    print("\n" + "="*35)

    # === Draw the Final Aircraft ===
    airplane_sol.draw(show=True)
