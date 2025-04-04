import os 
import pytest 
import sys 
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from mesh import Mesh 
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


@pytest.fixture
def mesh(domain):
    m = Mesh(domain, dx=0.1, dy=0.1)
    return m


def test_mesh_initialization(mesh):
    assert mesh.domain == mesh.domain 
    assert np.all(mesh.nodes == None)
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


def test_construct_nodes_has_correct_positions_2d(mesh):
    mesh.construct_nodes()
    assert mesh.nodes.shape == (6, 21)
    assert mesh.nodes[0, 0].position == (0.0, -0.25)
    assert mesh.nodes[1, 1].position == (.1, -.15)
    assert np.isclose(mesh.nodes[2, 2].position,(.2, -.05)).all()
    
def test_construct_nodes_has_correct_positions_1d(mesh):
    mesh.dy = None 
    mesh = Mesh(mesh.domain, dx=.1)
    mesh.construct_nodes()
    assert mesh.nodes.shape == (21,)
    assert mesh.nodes[0].position == (0.0, 0.0)
    assert mesh.nodes[1].position == (.1, 0.0)
    assert np.isclose(mesh.nodes[2].position,(.2, 0.0)).all()

def test_node_indexing_algebra_is_correct(mesh):
    mesh.construct_nodes()
    sum = np.array(mesh.nodes[0, 0].position) + np.array(mesh.nodes[1, 0].position)
    expected = np.array((0.0, -.4))
    assert np.isclose(sum, expected).all()
    p0 = np.array(mesh.nodes[0, 0].position)
    p1 = np.array(mesh.nodes[1, 0].position)
    delta = p1 - p0
    expected_delta = np.array([0.0, mesh.dy])  # your expected grid spacing in y

    assert np.allclose(delta, expected_delta)

    dx = mesh.dx
    p0 = np.array(mesh.nodes[2, 0].position)
    p1 = np.array(mesh.nodes[2, 1].position)
    delta_x = p1 - p0
    expected = np.array([dx, 0.0])
    assert np.allclose(delta_x, expected)

def test_mesh_shape_matches_expected(mesh):
    mesh.construct_nodes()
    expected_nx = int((mesh.domain.get_bounds()[0][1] - mesh.domain.get_bounds()[0][0]) / mesh.dx) + 1
    expected_ny = int((mesh.domain.get_bounds()[1][1] - mesh.domain.get_bounds()[1][0]) / mesh.dy) + 1
    assert mesh.nodes.shape == (expected_ny, expected_nx)


def test_all_nodes_inside_domain(mesh):
    mesh.construct_nodes()
    for row in mesh.nodes:
        for node in np.atleast_1d(row):
            if node is not None:
                x, y = node.position
                assert mesh.domain.is_inside(x, y)


def test_corner_nodes_positions(mesh):
    mesh.construct_nodes()
    ny, nx = mesh.nodes.shape
    corners = [
        mesh.nodes[0, 0],          # bottom-left
        mesh.nodes[0, nx-1],       # bottom-right
        mesh.nodes[ny-1, 0],       # top-left
        mesh.nodes[ny-1, nx-1]     # top-right
    ]
    for node in corners:
        if node is not None:
            x, y = node.position
            assert mesh.domain.is_inside(x, y)

def test_1d_mesh_y_is_zero(mesh):
    mesh.dy = None
    mesh = Mesh(mesh.domain, dx=.1)
    mesh.construct_nodes()
    for node in mesh.nodes:
        if node is not None:
            assert np.isclose(node.position[1], 0.0)


def test_uniform_spacing_x_and_y(mesh):
    mesh.construct_nodes()
    ny, nx = mesh.nodes.shape
    for j in range(ny):
        for i in range(1, nx):
            left = mesh.nodes[j, i-1]
            right = mesh.nodes[j, i]
            if left is not None and right is not None:
                dx = right.position[0] - left.position[0]
                assert np.isclose(dx, mesh.dx)

    for j in range(1, ny):
        for i in range(nx):
            below = mesh.nodes[j-1, i]
            above = mesh.nodes[j, i]
            if below is not None and above is not None:
                dy = above.position[1] - below.position[1]
                assert np.isclose(dy, mesh.dy)


