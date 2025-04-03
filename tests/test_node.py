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
    assert node.P is None

def test_define_coefficients(node):
    node.define_coefficients(De=1, Fe=1, Dw=2, Fw=3, Dn=4, Fn=5, Ds=6, Fs=7, Sc=8, Sp=9, dx=.1, dy=.1)
    assert node.aE == 1
    assert node.aW == 5
    assert node.aN == 4
    assert node.aS == 13
    assert isclose(node.b, 0.08, rel_tol=1e-9)
    assert isclose(node.aP, 1 + 5 + 4 + 13 - 0.09, rel_tol=1e-9)
    