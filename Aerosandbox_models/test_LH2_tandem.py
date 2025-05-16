import pytest
import aerosandbox as asb
import aerosandbox.numpy as np
from LH2_tandem import *
import math
def test_wing_weight():
    W_TO = 2177  # N
    n_ult = 3
    A = 36
    S = 260  # m²
    taper_ratio = 0.5
    t_c = 0.12
    V_H = 30  # m/s
    expected = 1426
    result = wing_weight(W_TO, n_ult, A, S, taper_ratio, t_c, V_H)

    assert math.isclose(result, expected, rel_tol=1e-2), f"Expected {expected}, got {result}"

def test_avionics_weight():
    assert avionics_weight(1000) == 30
    assert avionics_weight(0) == 0
    assert avionics_weight(2000) == 60

def test_landing_gear_weight():
    assert landing_gear_weight(1000) == 50
    assert landing_gear_weight(0) == 0
    assert landing_gear_weight(2000) == 100

def test_motor_weight():
    result = motor_weight(10000)
    expected = 9.81
    assert math.isclose(result, expected, rel_tol=1e-3)

    # Zero power
    assert motor_weight(0) == 0

    # Large power
    result = motor_weight(50_000)
    expected = 49.05
    assert math.isclose(result, expected, rel_tol=1e-3)
    
test_wing_weight()
test_avionics_weight()
test_landing_gear_weight()
test_motor_weight()