import os 
import sys 
import pytest 
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from control_surface import ControlSurface

@pytest.fixture
def vertical_control_surface():
    cs = ControlSurface()
    cs.u = 1.0
    return cs

@pytest.fixture
def horizontal_control_surface():
    cs = ControlSurface()
    cs.v = 1.0
    return cs

def test_control_surface_initialization():
    cs = ControlSurface()
    assert cs.u is None
    assert cs.v is None
    assert cs.orientation is None
    assert cs.area is None
    assert cs.position is None
    assert cs.b is None


def test_vertical_control_surface(vertical_control_surface):
    cs = vertical_control_surface
    assert cs.u == 1.0
    assert cs.v is None
    assert cs.orientation == 'vertical'
    assert cs.area is None
    assert cs.position is None
    assert cs.b is None


def test_horizontal_control_surface(horizontal_control_surface):
    cs = horizontal_control_surface
    assert cs.v == 1.0
    assert cs.u is None
    assert cs.orientation == 'horizontal'
    assert cs.area is None
    assert cs.position is None
    assert cs.b is None

def test_calculate_x_coefficients(vertical_control_surface):
    cs = vertical_control_surface
    cs.u_old = 0.001
    cs.calculate_x_coefficients(dx=0.0125, dy=0.0025, De=.08, Fe=1, Dw=0.08, Fw=1, Dn=0.4, Fn=.1, Ds=0.4, Fs=0.1, pe=.001, pw=.001, alphau=.5, is_2d=True)
    
    assert cs.u_old == 0.001
    assert cs.aE == 0.0002
    assert cs.aW == .0027
    assert np.isclose(cs.aN,0.005)
    assert np.isclose(cs.aS, 0.00625)
    assert np.isclose(cs.aP , 0.01415)
    assert np.isclose(cs.aP/.5, 0.0283)
    assert np.isclose(cs.b, 1.415e-5)

def test_calculate_y_coefficients(horizontal_control_surface):
    cs = horizontal_control_surface
    cs.v_old = 0.001
    cs.calculate_y_coefficients(dx=0.0125, dy=0.0025, De=.08, Fe=.79616, Dw=0.08, Fw=.815735, Dn=0.4, Fn=.1, Ds=0.4, Fs=0.1, ps=.001, pn=.001, alphav=.5, is_2d=True)
    assert np.isclose(cs.aE, 0.0002)
    assert np.isclose(cs.aW , .0022393375)
    assert np.isclose(cs.aN, 0.005)
    assert np.isclose(cs.aS, 0.00625)
    assert np.isclose(cs.aP , 0.0136404)
    assert np.isclose(cs.aP/.5, 0.027281)
    assert np.isclose(cs.b , 1.36404e-5)