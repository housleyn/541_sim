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


def test_construct_control_surfaces(mesh):
    mesh.construct_nodes()
    mesh.construct_control_surfaces()
    assert mesh.u_faces.shape == (6, 20)
    assert mesh.v_faces.shape == (5, 21)


def test_control_surface_positions_u(mesh):
    mesh.construct_nodes()
    mesh.construct_control_surfaces()
    assert np.isclose(mesh.u_faces[3,0].position, (.05, 0.05)).all()
    assert np.isclose(mesh.u_faces[3, 1].position, (.15, .05)).all()
    assert np.isclose(mesh.u_faces[1, 1].position, (0.15, -0.15)).all()

def test_control_surface_positions_v(mesh):
    mesh.construct_nodes()
    mesh.construct_control_surfaces()
    assert np.isclose(mesh.v_faces[0, 0].position, (0.0, -0.2)).all()
    assert np.isclose(mesh.v_faces[1, 1].position, (0.1, -0.1)).all()
    assert np.isclose(mesh.v_faces[3, 2].position, (0.2, 0.1)).all()


def test_that_indice_for_surface_and_node_are_compatible(mesh):
    mesh.construct_nodes()
    mesh.construct_control_surfaces()

    ny, nx = mesh.shape

    # u_faces: center between node[j, i] and node[j, i+1]
    for j in range(ny):
        for i in range(nx - 1):
            node = mesh.nodes[j, i]
            face = mesh.u_faces[j, i]
            if node is not None and face is not None:
                dx_half = 0.5 * mesh.dx
                delta_x = abs(face.position[0] - node.position[0])
                assert np.isclose(delta_x, dx_half)

    # v_faces: center between node[j, i] and node[j+1, i]
    for j in range(ny - 1):
        for i in range(nx):
            node = mesh.nodes[j, i]
            face = mesh.v_faces[j, i]
            if node is not None and face is not None:
                dy_half = 0.5 * mesh.dy
                delta_y = abs(face.position[1] - node.position[1])
                assert np.isclose(delta_y, dy_half)
def test_node_neighbor_positions_are_correct(mesh):
    mesh.construct_nodes()
    ny, nx = mesh.shape

    for j in range(ny):
        for i in range(nx - 1):
            n1 = mesh.nodes[j, i]
            n2 = mesh.nodes[j, i + 1]
            if n1 is not None and n2 is not None:
                dx_actual = n2.position[0] - n1.position[0]
                assert np.isclose(dx_actual, mesh.dx)

    for j in range(ny - 1):
        for i in range(nx):
            n1 = mesh.nodes[j, i]
            n2 = mesh.nodes[j + 1, i]
            if n1 is not None and n2 is not None:
                dy_actual = n2.position[1] - n1.position[1]
                assert np.isclose(dy_actual, mesh.dy)
def test_u_face_is_between_neighbor_nodes(mesh):
    mesh.construct_nodes()
    mesh.construct_control_surfaces()
    ny, nx = mesh.shape

    for j in range(ny):
        for i in range(nx - 1):
            n1 = mesh.nodes[j, i]
            n2 = mesh.nodes[j, i + 1]
            u_face = mesh.u_faces[j, i]
            if n1 is not None and n2 is not None and u_face is not None:
                x_avg = 0.5 * (n1.position[0] + n2.position[0])
                assert np.isclose(u_face.position[0], x_avg)

def test_valid_indices_are_physically_adjacent(mesh):
    mesh.construct_nodes()
    mesh.build_valid_index_lists()

    if len(mesh.shape) == 1:
        # 1D case
        for idx in mesh.interior_node_indices:
            if idx > 0 and idx - 1 in mesh.interior_node_indices:
                left = mesh.nodes[idx - 1]
                curr = mesh.nodes[idx]
                if left is not None and curr is not None:
                    dx = curr.position[0] - left.position[0]
                    assert np.isclose(dx, mesh.dx)
    else:
        # 2D case
        ny, nx = mesh.nodes.shape
        for j, i in mesh.interior_node_indices:
            curr = mesh.nodes[j, i]
            if curr is None:
                continue

            # Check East neighbor
            if i + 1 < nx and (j, i + 1) in mesh.interior_node_indices:
                east = mesh.nodes[j, i + 1]
                if east is not None:
                    dx = east.position[0] - curr.position[0]
                    assert np.isclose(dx, mesh.dx)

            # Check West neighbor
            if i - 1 >= 0 and (j, i - 1) in mesh.interior_node_indices:
                west = mesh.nodes[j, i - 1]
                if west is not None:
                    dx = curr.position[0] - west.position[0]
                    assert np.isclose(dx, mesh.dx)

            # Check North neighbor
            if j + 1 < ny and (j + 1, i) in mesh.interior_node_indices:
                north = mesh.nodes[j + 1, i]
                if north is not None:
                    dy = north.position[1] - curr.position[1]
                    assert np.isclose(dy, mesh.dy)

            # Check South neighbor
            if j - 1 >= 0 and (j - 1, i) in mesh.interior_node_indices:
                south = mesh.nodes[j - 1, i]
                if south is not None:
                    dy = curr.position[1] - south.position[1]
                    assert np.isclose(dy, mesh.dy)
