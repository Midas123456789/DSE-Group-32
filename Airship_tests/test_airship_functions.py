import sys
sys.path.append('./')

from airship_functions import *

def test_calculate_lift():
    """
    Test the calculate_lift function.
    """
    # Test case 1: Basic test case
    WH0 = 9950  # weight of airship in lbf
    assert calculate_lift(WH0) == WH0, f"Expected {WH0}, got {calculate_lift(WH0)}"

def test_calculate_dynamic_pressure_max():
    """
    Test the calculate_dynamic_pressure function.
    """
    # Test case 1: Basic test case
    rhoSL = 0.002377  # slug/ft³
    Vmax = 76  # ft/s
    expected_qmax = 6.86  # lbf/ft²
    assert round(calculate_dynamic_pressure_max(rhoSL, Vmax),2) == expected_qmax, f"Expected {expected_qmax}, got {calculate_dynamic_pressure_max(rhoSL, Vmax)}"


def test_calculate_lift_coefficient():
    """
    Test the calculate_lift_coefficient function.
    """
    # Test case 1: Basic test case
    qmax = 6.86  # dynamic pressure in lbf/ft²
    WH0 = 9950  # weight of airship in lbf
    Vol = 1000000  # volume of the airship in ft³
    expected_CL_maxpower = 0.145
    assert round(calculate_lift_coefficient(qmax, WH0, Vol),3) == expected_CL_maxpower, f"Expected {expected_CL_maxpower}, got {calculate_lift_coefficient(rhoSL, Vmax, WH0, Vol)}"


def test_calculate_drag_maxpower():
    """
    Test the calculate_drag_maxpower function.
    """
    # Test case 1: Basic test case
    CD0 = 0.02747  # zero-lift drag coefficient [dimensionless]
    K = 0.869  # induced drag factor [dimensionless]
    CL_maxpower = 0.145  # lift coefficient at maximum power [dimensionless]
    qmax = 6.86  # dynamic pressure in lbf/ft²
    Vol = 1000000  # volume of the airship in ft³
    expected_D_maxpower = 3138
    assert round(calculate_drag_maxpower(CD0, K, CL_maxpower, qmax, Vol),0) == expected_D_maxpower, f"Expected {expected_D_maxpower}, got {calculate_drag_maxpower(CD0, K, CL_maxpower, qmax, Vol)}"


def test_calculate_power_per_engine():
    """
    Test the calculate_power_per_engine function.
    """
    # Test case 1: Basic test case
    Vmax = 76  # maximum velocity of the airship in ft/s
    D_maxpower = 3138  # drag force at maximum power in lbf
    n_eng = 0.75  # engine efficiency
    NE = 2  # number of engines
    expected_power_per_engine = 289
    assert round(calculate_power_per_engine(Vmax, D_maxpower, n_eng, NE),0) == expected_power_per_engine, f"Expected {expected_power_per_engine}, got {calculate_power_per_engine(Vmax, D_maxpower, n_eng, NE)}"

def test_calculate_speed_coefficient():
    """
    Test the calculate_speed_coefficient function.
    """
    # Test case 1: Basic test case
    Vmax = 76  # maximum velocity of the airship in ft/s
    rhoSL = 0.002377  # slug/ft³
    P_hp = 289
    n = 20
    expected_speed_coefficient = 0.624
    assert round(calculate_speed_coefficient(rhoSL, Vmax, P_hp, n),3) == expected_speed_coefficient, f"Expected {expected_speed_coefficient}, got {calculate_speed_coefficient(Vmax, rhoSL)}"


def test_calculate_propeller_advance_ratio():
    """
    Test the calculate_propeller_advance_ratio function.
    """
    # Test case 1: Basic test case
    C_S = 0.624  # speed coefficient [dimensionless]
    expected_propeller_advance_ratio = 0.349  # expected advance ratio
    assert round(calculate_propeller_advance_ratio(C_S),3) == expected_propeller_advance_ratio, f"Expected {expected_propeller_advance_ratio}, got {calculate_propeller_advance_ratio(C_S)}"

def test_calculate_propeller_diameter():
    """
    Test the calculate_propeller_diameter function.
    """
    # Test case 1: Basic test case
    Vmax = 76  # maximum velocity of the airship in ft/s
    J = 0.349  # propeller advance ratio [dimensionless]
    n = 20  # rotations per second [/s]
    expected_propeller_diameter = 10.9  # expected propeller diameter in ft
    assert round(calculate_propeller_diameter(J, n, Vmax),1) == expected_propeller_diameter, f"Expected {expected_propeller_diameter}, got {calculate_propeller_diameter(J, n, Vmax)}"

def test_calculate_propeller_efficiency():
    """
    Test the calculate_propeller_efficiency function.
    """
    # Test case 1: Basic test case
    C_S = 0.624  # speed coefficient [dimensionless]
    expected_propeller_efficiency = 0.609  # expected propeller efficiency
    assert round(calculate_propeller_efficiency(C_S), 3) == expected_propeller_efficiency, f"Expected {expected_propeller_efficiency}, got {calculate_propeller_efficiency(C_S)}"

def test_calculate_internal_presssure():
    """
    Test the calculate_internal_pressure function.
    """
    # Test case 1: Basic test case
    qmax = 6.86  # dynamic pressure in lbf/ft²
    height = 78.2 # height of the airship in ft
    expected_internal_pressure = 0.092  # expected internal pressure in psi
    assert round(calculate_internal_pressure(qmax, height), 3) == expected_internal_pressure, f"Expected {expected_internal_pressure}, got {calculate_internal_pressure(qmax, height)}"

def test_calculate_hull_fabric_load():
    """
    Test the hull_fabirc_load function.
    """
    # Test case 1: Basic test case
    FS = 4 # factor of safety [dimensionaless]
    P_int = 0.092  # internal pressure psi
    height = 78.2 # height of the airship in ft
    expected_hull_fabric_load = 173  # expected internal pressure in psi
    assert round(calculate_hull_fabric_load(FS, P_int, height), 0) == expected_hull_fabric_load, f"Expected {expected_hull_fabric_load}, got {calculate_hull_fabric_load(FS, P_int, height)}"


def test_calculate_hull_fabric_density():
    """
    Test the calculate_hull_fabric_density function.
    """
    # Test case 1: Basic test case
    material = "Polyester (weave)"  # material type
    q_hull = 173  # hull fabric load in lb/in
    expected_hull_fabric_density = 9.80  # expected hull fabric density in oz/yd², differs from the book (states 9.75)
    assert round(calculate_hull_fabric_density(material, q_hull),
                 2) == expected_hull_fabric_density, f"Expected {expected_hull_fabric_density}, got {calculate_hull_fabric_density(material, q_hull)}"


def test_calculate_weight_envelope():
    """
    Test the calculate_weight_envelope function.
    """
    # Test case 1: Basic test case
    w_hull = 9.75  # hull fabric density in oz/yd²
    S_wet = 61631  # wetted area in ft²
    expected_weight = 6309 # expected weight of envelope in lb, differs from the book (states 6311)

    assert round(calculate_weight_envelope(w_hull, S_wet),
                 0) == expected_weight, f"Expected {expected_weight}, got {calculate_weight_envelope(w_hull, S_wet)}"


def test_calculate_septum_density():
    """
    Test the calculate_septum_density function.
    """
    # Test case 1: Basic test case
    material = "Polyester (weave)"  # material type
    q_hull = 173  # septum fabric load in lb/in
    expected_septum_density = 13.72  # expected septum density in oz/yd²
    assert round(calculate_septum_density(material, q_hull),
                 2) == expected_septum_density, f"Expected {expected_septum_density}, got {calculate_septum_density(material, q_hull)}"
