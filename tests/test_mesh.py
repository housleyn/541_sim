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
    m = Mesh(domain, nx=20, ny=2)
    return m


def test_interior_node_indices(mesh):
    # Construct the mesh
    mesh.construct_nodes()
    
    # For 1D mesh, interior_node_indices should be a list of single integers
    if len(mesh.shape) == 1:
        assert all(isinstance(idx, int) for idx in mesh.interior_node_indices), \
            "In a 1D mesh, interior_node_indices should be a list of integers."
    
    # For 2D mesh, interior_node_indices should be a list of tuples (j, i)
    elif len(mesh.shape) == 2:
        assert all(isinstance(idx, tuple) and len(idx) == 2 for idx in mesh.interior_node_indices), \
            "In a 2D mesh, interior_node_indices should be a list of tuples (j, i)."
        
        # Additionally, check if the indices (j, i) are valid and within the mesh shape
        for j, i in mesh.interior_node_indices:
            assert 0 <= j < mesh.shape[0], f"Invalid row index {j}."
            assert 0 <= i < mesh.shape[1], f"Invalid column index {i}."
    
    # Print out the result for debugging (optional)
    print("Interior node indices:", mesh.interior_node_indices)


def test_mesh_initialization(mesh):
    assert mesh.domain == mesh.domain
    assert np.all(mesh.nodes == None)
    assert np.allclose(mesh.dx, 0.1)

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
    assert mesh.shape == (2, 20)
    mesh.dy = None 
    mesh = Mesh(mesh.domain, nx=20)
    assert mesh.shape == (20,)


def test_construct_nodes_has_correct_positions_2d(mesh):
    mesh.construct_nodes()

    assert mesh.nodes.shape == (2, 20)

    node00 = mesh.nodes[0, 0]
    node01 = mesh.nodes[0, 1]
    node10 = mesh.nodes[1, 0]

    assert node00 is not None
    assert np.isclose(node00.position[0], 0.0)

    if node01 is not None:
        assert np.isclose(node01.position[0], mesh.dx[1])

    if node10 is not None:
        assert np.isclose(node10.position[0], 0.0)



    
def test_construct_nodes_has_correct_positions_1d(mesh):
    mesh.dy = None 
    mesh = Mesh(mesh.domain, nx=20)
    mesh.construct_nodes()
    assert mesh.nodes.shape == (20,)
    assert mesh.nodes[0].position == (0.0, 0.0)
    assert mesh.nodes[1].position == (.1, 0.0)
    assert np.isclose(mesh.nodes[2].position,(.2, 0.0)).all()

def test_node_indexing_algebra_is_correct(mesh):
    mesh.construct_nodes()

    node0 = mesh.nodes[0, 0]
    node1 = mesh.nodes[1, 0]

    assert node0 is not None and node1 is not None

    dx_actual = node1.position[0] - node0.position[0]
    dy_actual = node1.position[1] - node0.position[1]

    assert np.isclose(dx_actual, 0.0)  # Same x-column
    assert np.isclose(dy_actual, mesh.dy[1])  # Vertical spacing based on mesh.dy


def test_mesh_shape_matches_expected(mesh):
    mesh.construct_nodes()
    assert mesh.nodes.shape == (mesh.ny, mesh.nx)



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
    mesh = Mesh(mesh.domain, nx=20)
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
                assert np.isclose(dx, mesh.dx[i-1])

    for j in range(1, ny):
        for i in range(nx):
            below = mesh.nodes[j-1, i]
            above = mesh.nodes[j, i]
            if below is not None and above is not None:
                dy = above.position[1] - below.position[1]
                assert np.isclose(dy, mesh.dy[j-1])



def test_construct_control_surfaces(mesh):
    mesh.construct_nodes()
    mesh.construct_control_surfaces()
    assert mesh.u_faces.shape == (2, 19)
    assert mesh.v_faces.shape == (1, 20)


def test_control_surface_positions_u(mesh):
    mesh.construct_nodes()
    mesh.construct_control_surfaces()

    # Ensure there's a control surface at (1, 0) in the middle row
    face = mesh.u_faces[1, 0]  # Check in the middle row instead of the bottom row

    # Check if the face is not None
    assert face is not None, "u_face at (1, 0) should not be None"

    # Check that mesh.nodes[0, 0] and mesh.nodes[1, 0] are not None
    assert mesh.nodes[1, 0] is not None, "Node at (1, 0) should not be None"
    assert mesh.nodes[1, 1] is not None, "Node at (1, 1) should not be None"

    # Now check if the position is correct for the u face at (1, 0)
    expected_x = 0.5 * (mesh.nodes[1, 0].position[0] + mesh.nodes[1, 1].position[0])  # Midpoint of x positions in row 1
    assert np.isclose(face.position[0], expected_x), f"Expected x position: {expected_x}, but got: {face.position[0]}"

    # Optionally check the y-position of the face if needed
    expected_y = 0.5 * (mesh.nodes[1, 0].position[1] + mesh.nodes[1, 1].position[1])  # Midpoint of y positions
    assert np.isclose(face.position[1], expected_y), f"Expected y position: {expected_y}, but got: {face.position[1]}"




def test_control_surface_positions_v(mesh):
    mesh.construct_nodes()
    mesh.construct_control_surfaces()

    # Check the position of the first v-face
    face = mesh.v_faces[0, 0]

    # Ensure the face is not None
    assert face is not None, "v_face at (0, 0) should not be None"

    # Calculate the expected position based on the first two nodes' y-values
    expected_y = 0.5 * (mesh.nodes[0, 0].position[1] + mesh.nodes[1, 0].position[1])

    # Now, check if the calculated y position matches the expected value
    assert np.isclose(face.position[1], expected_y), f"Expected y position: {expected_y}, but got: {face.position[1]}"

    # Optionally check the x-position if needed (it should be consistent with the x of nodes)
    expected_x = mesh.nodes[0, 0].position[0]
    assert np.isclose(face.position[0], expected_x), f"Expected x position: {expected_x}, but got: {face.position[0]}"



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
                dx_half = 0.5 * mesh.dx[i]  # Use the specific dx value for this index i
                delta_x = abs(face.position[0] - node.position[0])  # Ensure scalar comparison

                assert np.isclose(delta_x, dx_half), f"Expected delta_x: {dx_half}, but got: {delta_x}"

def test_node_neighbor_positions_are_correct(mesh):
    mesh.construct_nodes()
    ny, nx = mesh.shape
    
    for j in range(ny):
        for i in range(nx - 1):
            n1 = mesh.nodes[j, i]
            n2 = mesh.nodes[j, i + 1]
            if n1 is not None and n2 is not None:
                dx_actual = n2.position[0] - n1.position[0]
                dx_expected = mesh.dx[i]  # Use the correct dx value for this cell

                assert np.isclose(dx_actual, dx_expected), f"Expected dx: {dx_expected}, but got: {dx_actual}"

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
                    assert np.isclose(dx, mesh.dx[idx])  # Index into mesh.dx to get the correct dx for the node
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
                    assert np.isclose(dx, mesh.dx[i])  # Access correct dx for the node

            # Check West neighbor
            if i - 1 >= 0 and (j, i - 1) in mesh.interior_node_indices:
                west = mesh.nodes[j, i - 1]
                if west is not None:
                    dx = curr.position[0] - west.position[0]
                    assert np.isclose(dx, mesh.dx[i])  # Access correct dx for the node

            # Check North neighbor
            if j + 1 < ny and (j + 1, i) in mesh.interior_node_indices:
                north = mesh.nodes[j + 1, i]
                if north is not None:
                    dy = north.position[1] - curr.position[1]
                    assert np.isclose(dy, mesh.dy[j])  # Access correct dy for the node

def test_mesh_has_nonempty_interior_indices(mesh):
    # Ensure that mesh has been constructed and valid index lists built
    mesh.construct_nodes()
    mesh.construct_control_surfaces()
    mesh.build_valid_index_lists()

    # Check that interior_node_indices is not empty
    assert len(mesh.interior_node_indices) > 0, "interior_node_indices is empty but should not be."
