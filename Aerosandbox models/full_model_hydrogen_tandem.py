import aerosandbox as asb
import aerosandbox.numpy as np
from mass_wing import Mass_wing

def get_LH_conventional(N_cords = 5, wing_airfoil = asb.Airfoil("sd7037"), altitude = 18000, mission_days = 10, W_payload = 1000, P_payload = 10e3):
    g = 9.81
    atm = asb.Atmosphere(altitude=altitude)
    opti = asb.Opti()
    Tail_drag_factor = 1.03
    second_in_day = 86400
    mission_seconds = mission_days * second_in_day

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
    wing1 = Mass_wing(
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

    wing2 = Mass_wing(
        symmetric=True,
        xsecs=[
            asb.WingXSec(
                xyz_le=[-0.25 * cords[i], y_sections[i], 0],
                chord=cords[i],
                airfoil=wing_airfoil
            )
            for i in range(N_cords)
        ]
    ).translate([7, 0, 0])

    wings = [wing1] + [wing2]
    airplane = asb.Airplane(wings=wings)

    # === Aerodynamics === #
    op = asb.OperatingPoint(velocity=V, alpha=alpha, atmosphere=atm)
    aero = asb.AeroBuildup(airplane=airplane, op_point=op).run()
    
    total_area = wing1.area() + wing2.area()
    
    # === Lift & Drag Forces === #
    L = 0.5 * atm.density() * V**2 * aero['CL'] * total_area
    D = 0.5 * atm.density() * V**2 * aero['CD'] * total_area * Tail_drag_factor

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
    W_skin = 2 * total_area * skin_density * g

    # === Wing Weight === #
    spar_mass, t_list = wing1.spar_mass(L/2, b)
    spar_mass_total = 2 * spar_mass
    W_spar = spar_mass_total * g   
    W_wing = W_spar + W_skin
    M_fuel_cell = E_LH / mission_seconds / power_density_PEMFCs
    W_fuel_cell = M_fuel_cell * g
    # === Miss Weight === #
    W_misc = W_wing

    # === Total Weight === #
    W_total = W_payload + W_wing + W_LH + W_misc + W_tank + W_fuel_cell
    opti.subject_to(L >= W_total)
    
    opti.subject_to([
    cords > 0,
    np.diff(cords) <= 0,  # Taper
    total_area < 2000
    ])

    # === Objective === #
    opti.minimize(P_required)
    # === Solve === #
    sol = opti.solve(verbose=False)
    CD = aero["CD"] * Tail_drag_factor
    CL = aero["CL"]
    return {
    "sol": sol,
    "cords": cords,
    "b": b,
    "V": V,
    "alpha": alpha,
    "wing_area": total_area,
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
    "M_LH": M_LH,
    "M_fuel_cell": M_fuel_cell,
    "M_tank": M_tank,
    "LH_volume": LH_volume,
    "Airfoil": wing_airfoil,
    "Altitude": altitude,
    "Atm": atm,
    "W_payload": W_payload,
    "P_payload": P_payload,
    
    } 


if __name__ == "__main__":
    results = get_LH_conventional()
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
    M_LH_sol = sol.value(results["M_LH"])
    M_fuel_cell_sol = sol.value(results["M_fuel_cell"])
    M_tank_sol = sol.value(results["M_tank"])
    LH_volume_sol = sol.value(results["LH_volume"])
    altitude = sol.value(results["Altitude"])
    atm = sol.value(results["Atm"])
    wing_airfoil = sol.value(results["Airfoil"])
    W_payload = sol.value(results["W_payload"])
    P_payload = sol.value(results["P_payload"])
    W_LH = M_LH_sol * 9.81
    Total_weight_sol = W_total_sol

    # === Rebuild Solved Wing ===
    y_sections_sol = np.linspace(0, b_sol / 2, len(cords_sol))
    wing_sol1 = Mass_wing(
        symmetric=True,
        xsecs=[
            asb.WingXSec(
                xyz_le=[-0.25 * cords_sol[i], y_sections_sol[i], 0],
                chord=cords_sol[i],
                airfoil=wing_airfoil
            )
            for i in range(len(cords_sol))
        ]
    ).translate([0, 0, 0])

    wing_sol2 = Mass_wing(
        symmetric=True,
        xsecs=[
            asb.WingXSec(
                xyz_le=[-0.25 * cords_sol[i], y_sections_sol[i], 0],
                chord=cords_sol[i],
                airfoil=wing_airfoil
            )
            for i in range(len(cords_sol))
        ]
    ).translate([9, 0, 0])

    airplane_sol = asb.Airplane(wings=[wing_sol1] + [wing_sol2])

    # === Output Results ===
    print("\n" + "="*35)
    print("        OPTIMIZATION RESULTS")
    print("="*35 + "\n")

    print("=== Geometry ===")
    print(f"  Span (b):                 {b_sol:.2f} m")
    print(f"  Chord distribution (m):   {cords_sol}")
    print(f"  Wing area:                {area_sol:.2f} m²")
    print(f"  Volume LH:                {LH_volume_sol:.2f} m³\n")

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
    print(f"  LH weight:                {W_LH:.2f} N")
    print(f"  Miscellaneous weight:     {W_misc_sol:.2f} N")
    print(f"  Total weight:             {Total_weight_sol:.2f} N\n")

    print("=== Mass ===")
    print(f"  Payload mass:             {W_payload / 9.81:.2f} kg")
    print(f"  Spar mass:                {W_spar_sol / 9.81:.2f} kg")
    print(f"  Skin mass:                {W_skin_sol / 9.81:.2f} kg")
    print(f"  Miscellaneous mass:       {W_misc_sol / 9.81:.2f} kg")
    print(f"  LH mass:                  {M_LH_sol:.2f} kg")
    print(f"  Fuel cell mass:           {M_fuel_cell_sol:.2f} kg")
    print(f"  Tank mass:                {M_tank_sol:.2f} kg")
    print(f"  Total mass:               {Total_weight_sol / 9.81:.2f} kg\n")

    print("=== Operations ===")
    print(f"  Yearly LH consumption:    {M_LH_sol * 52:.2f} kg")
    eu_per_kg_LH = 10
    print(f"  Yearly Fuel Cost:         {M_LH_sol * 52 * eu_per_kg_LH:.2f} Euro\n")

    print("=== Structural ===")
    print(f"  Spar thicknesses (t_list):\n    {t_list_sol}")

    print("\n" + "="*35)

    # === Draw the Final Aircraft ===
    airplane_sol.draw(show=True)
