import os 
import pytest 
import sys 
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from mesh import Mesh 
from domain import Domain

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
    m = Mesh(domain, dx=0.1, dy=0.1)
    return m


def test_mesh_initialization(mesh):
    assert mesh.domain == mesh.domain 
    assert mesh.nodes == [] 
    assert mesh.dx == .1 
    assert mesh.dy == .1
    mesh.dy = None
    assert mesh.dy == None
    mesh.dy = .1

def test_mesh_reads_domain(mesh):
    assert mesh.domain.lower_boundary_func(0, 0) == 0.25
    assert mesh.domain.upper_boundary_func(0, 0) == -0.25
    assert mesh.domain.left_boundary_func(0) == 0
    assert mesh.domain.right_boundary_func(0) == 2
    assert mesh.domain.is_inside(1, 0) == True
    assert mesh.domain.is_inside(0, 0) == True
    assert mesh.domain.is_inside(2, 0) == True
    assert mesh.domain.is_inside(0, -1) == False
    assert mesh.domain.is_inside(2, -1) == False
    assert mesh.domain.is_inside(3, 0) == False
    assert mesh.domain.is_inside(-1, 0) == False
    assert mesh.domain.is_inside(0, .2) == True #on boundary
    assert mesh.domain.is_inside(1, .15) == True #on boundary

def test_shape_constructs_correctly_for_1D_and_2D(mesh):
    assert mesh.shape == (6, 21)
    mesh.dy = None 
    mesh = Mesh(mesh.domain, dx=.1)
    assert mesh.shape == (21,)




