import os 
import pytest 
import sys 
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from node import Node
from math import isclose

@pytest.fixture 
def node():
    return Node()

def test_node_class_initialized_correctly(node):
    assert node.aW is None
    assert node.aE is None
    assert node.aP is None 
    assert node.aN is None 
    assert node.aS is None 
    assert node.b is None 
    assert node.p is None

def test_define_pressure_correction_coefficients(node):
    rho = 1000
    dx = 0.0125
    dy = 0.0025
    aP_u_E = .01415
    aP_u_W = .01415
    u_star_W = .01415
    u_star_E = .01415
    v_star_S = 1.0
    v_star_N = 2.0

    aP_v_N = 0.005
    aP_v_S = 0.00625
    alphau = .5
    alphav = .5

    node.define_pressure_correction_coefficients(rho, alphau, alphav,dx, dy, aP_u_W, aP_u_E, aP_v_N, aP_v_S,
                                                  u_star_W, u_star_E, v_star_S, v_star_N)

    assert isclose(node.aE, 0.2208480565371025)
    