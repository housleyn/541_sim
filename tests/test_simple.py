import os 
import sys
import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from mesh import Mesh 
from domain import Domain 
from material import Material
from simple import SIMPLE
from boundaries import Boundary

@pytest.fixture
def domain():
    d = Domain()
    d.define_lower_boundary(lambda x, y: y - .1 * x + 0.25)
    d.define_upper_boundary(lambda x, y: y + .1 * x - 0.25)
    d.define_left_boundary(lambda y, x=None:0)
    d.define_right_boundary(lambda y, x=None:2)
    return d
@pytest.fixture
def mesh(domain):
    mesh = Mesh(domain, dx=0.1, dy=0.1)
    return mesh
@pytest.fixture
def boundary(mesh):
    bc = Boundary(mesh)
    return bc
@pytest.fixture
def material():
    m = Material()
    m.rho = 1.0
    m.mu = 0.0
    return m
@pytest.fixture
def simple(mesh, boundary, material):
    s = SIMPLE(mesh, boundary, material)
    return s

def test_all_nodes_have_information_set(simple):
    for i in range(simple.mesh.nodes.shape[0]):
        for j in range(simple.mesh.nodes.shape[1]):
            node = simple.mesh.nodes[i, j]
            if node is not None:
                assert node.aP is not None
                assert node.aE is not None
                assert node.aW is not None
                assert node.aN is not None
                assert node.aS is not None
                assert node.b is not None

def test_all_u_faces_have_information_set(simple):
    for i in range(simple.mesh.u_faces.shape[0]):
        for j in range(simple.mesh.u_faces.shape[1]):
            face = simple.mesh.u_faces[i, j]
            if face is not None:
                assert face.u is not None 

def test_all_v_faces_have_information_set(simple):
    for i in range(simple.mesh.v_faces.shape[0]):
        for j in range(simple.mesh.v_faces.shape[1]):
            face = simple.mesh.v_faces[i, j]
            if face is not None:
                assert face.v is not None