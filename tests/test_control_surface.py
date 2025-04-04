import os 
import sys 
import pytest 

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
