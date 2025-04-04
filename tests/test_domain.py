import os 
import pytest 
import sys 
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from domain import Domain
import numpy as np 


@pytest.fixture
def domain():
    d = Domain()
    d.define_lower_boundary(lambda x, y: y - .1 * x + 0.25)
    d.define_upper_boundary(lambda x, y: y + .1 * x - 0.25)
    d.define_left_boundary(lambda y, x=None:0)
    d.define_right_boundary(lambda y, x=None:2)
    return d


def test_lower_boundary_is_correct(domain):
    x_vals = np.linspace(0,2,100)
    expected = -0.1 * x_vals + .25 
    result = np.array([domain.lower_boundary_func(x, 0) for x in x_vals])
    assert np.allclose(result, expected, rtol=1e-5)

def test_upper_boundary_is_correct(domain):
    x_vals = np.linspace(0,2,100)
    expected = 0.1 * x_vals - .25
    result = np.array([domain.upper_boundary_func(x, 0) for x in x_vals])
    assert np.allclose(result, expected, rtol=1e-5)

def test_left_boundary_is_correct(domain):
    y_vals = np.linspace(-.5, .5, 100)
    expected = 0
    result = np.array([domain.left_boundary_func(y) for y in y_vals])
    assert np.allclose(result, expected, rtol=1e-5)

def test_right_boundary_is_correct(domain):
    y_vals = np.linspace(-.5, .5, 100)
    expected = 2
    result = np.array([domain.right_boundary_func(y) for y in y_vals])
    assert np.allclose(result, expected, rtol=1e-5)