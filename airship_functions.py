# Exercise/Line 19
def calculate_lift(WH0):
    """
    Calculate the lift force of the airship based on the weight of the airship and the buoyancy force.
    Inputs: WHO (float): weight of airship in lbf
    Outputs: lift (float): lift force in lbf
    """
    L_a = WH0
    return L_a

# Exercise/Line 20
def calculate_dynamic_pressure_max(rhoSL, Vmax):
    """
    Calculate the maximum dynamic pressure of the airship based on the density of air and maximum velocity.
    Inputs: rhoSL (float): density of air in slug/ft³
            Vmax (float): maximum velocity of the airship in ft/s
    Outputs: qmax (float): maximum dynamic pressure in lbf/ft²
    """
    # rhoSL = 0.002377  # slug/ft³

    qmax = 0.5 * rhoSL * Vmax**2  # dynamic pressure in lbf/ft²
    return qmax


def calculate_lift_coefficient(qmax, WH0, Vol):
    """
    Calculate the lift coefficient of the airship based on the weight of the airship and the buoyancy force.
    Inputs: qmax (float): maximum dynamic pressure in lbf/ft²
            WH0 (float): weight of airship in lbf
            Vol (float): volume of the airship in ft³
    Outputs: CL_maxpower (float): lift coefficient [dimensionless]
    """

    CL_maxpower = WH0/(qmax * Vol**(2/3)) # lift coefficient at maximum power

    return CL_maxpower


# Exercise/Line 21
def calculate_drag_maxpower(CD0, K, CL_maxpower, qmax, Vol):
    """
    Calculate the drag force of the airship at maximum power based on the drag coefficient, lift coefficient, and dynamic pressure.
    Inputs: CD0 (float): zero-lift drag coefficient [dimensionless]
            K (float): induced drag factor [dimensionless]
            CL_maxpower (float): lift coefficient at maximum power [dimensionless]
            qmax (float): dynamic pressure in lbf/ft²
            Vol (float): volume of the airship in ft³
    Outputs: D_maxpower (float): drag force at maximum power in lbf

    """

    D_maxpower = (CD0 + K * CL_maxpower**2) * qmax * Vol**(2/3)  # drag force at maximum power in lbf

    return D_maxpower

# Exercise/Line 22
def calculate_power_per_engine(Vmax, D_maxpower, n_eng, NE):
    """
    Calculate the power required per engine at maximum velocity.
    Inputs: Vmax (float): maximum velocity of the airship in ft/s
            D_maxpower (float): drag force at maximum power in lbf
            n_eng (float): engine efficiency [dimensionless]
            NE (int): number of engines [#]
    Outputs: P_per_engine (float): power required per engine in hp

    """

    # 550 = 550 ft-lbf/s = 1 hp
    P_hp_per_engine = ( (D_maxpower * Vmax) / (n_eng * NE)) / 550  # power required per engine in hp

    return P_hp_per_engine

# NOTE Skipped 23

# Exercise/Line 24
def calculate_speed_coefficient(rho, Vmax, P_hp, n):
    """
    Calculate the speed coefficient of the airship based on the density of air, maximum velocity, power, and number of engines.
    Inputs: rho (float): density of air in slug/ft³
            Vmax (float): maximum velocity of the airship in ft/s
            P (float): power per engine in hp
            n (float): rotations per second [/s]

    Outputs: C_speed (float): speed coefficient [dimensionless]

    """
    P = P_hp * 550  # Convert power from hp to ft-lbf/s

    C_S = ((rho * Vmax**5) / (P * n**2))**(1/5)  # speed coefficient [dimensionless]

    return C_S

def calculate_propeller_advance_ratio(C_S):
    """
    Calculate the propeller advance ratio based on the speed coefficient.
    Inputs: C_S (float): speed coefficient [dimensionless]
    Outputs: J (float): propeller advance ratio [dimensionless]

    """
    # NOTE Relation Clark Y airfoil, figure 5.9a
    J = 0.156 * C_S**2 + 0.241 * C_S + 0.138  # propeller advance ratio [dimensionless]

    return J


def calculate_propeller_diameter(J, n, Vmax):
    """
    Calculate the propeller diameter based on the advance ratio and speed coefficient.
    Inputs: J (float): propeller advance ratio [dimensionless]
            n (float): rotations per second [/s]
            Vmax (float): maximum velocity of the airship in ft/s
    Outputs: D_prop (float): propeller diameter in ft

    """
    D_prop = Vmax/(J*n)  # propeller diameter in ft

    return D_prop

# Exercise/Line 25

def calculate_propeller_efficiency(C_S):
    """
    Calculate the propeller efficiency based on the speed coefficient.
    Inputs: C_S (float): speed coefficient [dimensionless]
    Outputs: eta_prop (float): propeller efficiency [dimensionless]

    """
    # NOTE Relation Clark Y airfoil, figure 5.9a
    eta_prop = 0.139 * C_S**3 -0.749 * C_S**2 +1.37*C_S + 0.0115  # propeller efficiency [dimensionless]

    return eta_prop


# Exercise/Line 26
def calculate_internal_pressure(qmax, height):
    """
    Calculate the internal pressure of the airship based on the dynamic pressure and height.
    Inputs: qmax (float): maximum dynamic pressure in lbf/ft²
            height (float): height of the airship in ft
    Outputs: P_int (float): internal pressure in psi

    """

    P_int = 1.2*qmax + 0.0635*height  # internal pressure in lbf/ft²
    P_int = P_int / 144  # Convert from lbf/ft² to psi
    # NOTE 1 psi = 144 lbf/ft²

    return P_int

# Exercise/Line 27
def calculate_hull_fabric_load(FS, P_int, height):
    """
    Calculate the hull fabric load based on the safety factor, internal pressure, and radius of the hull.
    Inputs: FS (float): safety factor [dimensionless]
            P_int (float): internal pressure in psi
            height (float): height of the hull in ft
    Outputs: q_hull (float): hull fabric load in lb/in
            q for shear flow in lb/in

    """
    height_inches = height * 12  # Convert height from ft to inches
    # NOTE 1 ft = 12 inches
    q_hull = FS * P_int * height_inches  # hull fabric load in psi

    return q_hull



