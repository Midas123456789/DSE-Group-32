import aerosandbox as asb
import aerosandbox.numpy as np
from mass_wing import Mass_wing

def get_LH_flying_wing(N_cords = 5, wing_airfoil = asb.Airfoil("sd7037"), altitude = 18000, mission_days = 10, W_payload = 5e3, P_payload = 10e3):
    g = 9.81
    atm = asb.Atmosphere(altitude=altitude)
    opti = asb.Opti()
    Tail_drag_factor = 1.02
    second_in_day = 86400
    mission_days = 10
    mission_seconds = mission_days * second_in_day
    density_LH = 70.85 #kg/m3
    energy_density_LH = 142 * 10 ** 6 #J/kg
    power_density_PEMFCs = 40000 #W/kg
    PEMFCs_eff = 0.55
    motor_eff = 0.97 #Fundamentals of Aircraft and Airship Design
    propeller_eff = 0.85 #Fundamentals of Aircraft and Airship Design
    conversion_eff = PEMFCs_eff
    propulsive_eff = conversion_eff * motor_eff * propeller_eff
    fuel_tank_fuel_mass_fraction = 0.34  # from Brewer, Hydrogen Aircraft Technology pg. 29

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
    
    opti.subject_to([
    cords > 0,
    np.diff(cords) <= 0,  # Taper
    wing.area() < 2000
    ])
    
    wing_area = wing.area()
    CL = aero['CL']
    CD = aero['CD'] * Tail_drag_factor
    
    # === Objective === #
    opti.minimize(P_required)
    # === Solve === #
    try:
        sol = opti.solve(verbose=False)
        return {
        "cords": sol.value(cords),
        "b": sol.value(b),
        "V": sol.value(V),
        "alpha": sol.value(alpha),
        "wing_area": sol.value(wing.area()),
        "CL": sol.value(CL),
        "CD": sol.value(CD),
        "L": sol.value(L),
        "D": sol.value(D),
        "P_required": sol.value(P_required),
        "W_spar": sol.value(W_spar),
        "t_list": sol.value(t_list),
        "W_total": sol.value(W_total),
        "W_misc": sol.value(W_misc),
        "W_skin": sol.value(W_skin),
        "M_LH": sol.value(M_LH),
        "M_fuel_cell": sol.value(M_fuel_cell),
        "M_tank": sol.value(M_tank),
        "LH_volume": sol.value(LH_volume),
        "Airfoil": wing_airfoil,  # Only keep this if wing_airfoil is not symbolic
        "Altitude": altitude,      # Assumed numeric input
        "Atm": atm                 # Only if atm is not symbolic
        }
    except Exception as e:
        print(f"No solution found. Error: {e}")
        return

if __name__ == "__main__":
    # === Unpack Results from Dictionary ===
    results = get_LH_flying_wing()

    cords_sol = results["cords"]
    b_sol = results["b"]
    V_sol = results["V"]
    alpha_sol = results["alpha"]
    area_sol = results["wing_area"]
    CL_sol = results["CL"]
    CD_sol = results["CD"]
    L_sol = results["L"]
    D_sol = results["D"]
    P_required_sol = results["P_required"]
    W_spar_sol = results["W_spar"]
    t_list_sol = results["t_list"]
    W_total_sol = results["W_total"]
    W_misc_sol = results["W_misc"]
    W_skin_sol = results["W_skin"]
    M_LH_sol = results["M_LH"]
    M_fuel_cell_sol = results["M_fuel_cell"]
    M_tank_sol = results["M_tank"]
    LH_volume_sol = results["LH_volume"]
    altitude = results["Altitude"]
    atm = results["Atm"]
    wing_airfoil = results["Airfoil"]

    
    W_LH = M_LH_sol * 9.81
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
    #print(f"  Payload power:            {P_payload:.2f} W")
    print(f"  Power required (P):       {P_required_sol:.2f} W\n")

    print("=== Weights ===")
   # print(f"  Payload weight:           {W_payload:.2f} N")
    print(f"  Spar weight:              {W_spar_sol:.2f} N")
    print(f"  Skin weight:              {W_skin_sol:.2f} N")
    print(f"  LH weight:                {W_LH:.2f} N")
    print(f"  Miscellaneous weight:     {W_misc_sol:.2f} N")
    print(f"  Total weight:             {Total_weight_sol:.2f} N\n")

    print("=== Mass ===")
    #print(f"  Payload mass:             {W_payload / 9.81:.2f} kg")
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
