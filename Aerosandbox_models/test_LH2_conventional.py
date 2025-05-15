import pytest
import aerosandbox as asb
import aerosandbox.numpy as np
from LH2_conventional import Conventional_LH


@pytest.fixture
def lh_model():
    model = Conventional_LH(N_cords=3)  # Keep small for test speed
    model._variables_setup()
    model._setup_op()
    return model


def test_variable_setup(lh_model):
    assert lh_model.cords.shape == (3,1)
    assert lh_model.b.shape == (1,1)
    assert lh_model.V.shape == ()
    assert lh_model.alpha.shape == ()
    assert lh_model.op is not None


def test_geometry_setup(lh_model):
    lh_model._plane_geom_setup()
    assert len(lh_model.wing.xsecs) == 3
    assert isinstance(lh_model.airplane, asb.Airplane)


def test_aero_setup(lh_model):
    lh_model._plane_geom_setup()
    lh_model._setup_aero()
    assert "CL" in lh_model.aero
    assert "CD" in lh_model.aero
    assert lh_model.L is not None
    assert lh_model.D is not None


def test_compute_weights(lh_model):
    lh_model._plane_geom_setup()
    lh_model._setup_aero()
    lh_model.compute_weights()
    assert lh_model.W_total > 0
    assert lh_model.LH_volume > 0
    assert lh_model.W_wing > 0
    assert lh_model.W_LH > 0
    assert lh_model.W_misc > 0
    assert lh_model.W_spar > 0
    assert lh_model.W_skin > 0


def test_constraints_setup(lh_model):
    lh_model._plane_geom_setup()
    lh_model._setup_aero()
    lh_model.compute_weights()
    lh_model._setup_constraints()
    # No assert needed — if constraints break bounds, AeroSandbox will error during solve


def test_solver_runs():
    model = Conventional_LH(N_cords=3)
    sol = model.solve()
    assert "b" in sol
    assert sol["W_total"] > 0
    assert sol["wing_area"] > 0
    assert sol["P_required"] > 0
