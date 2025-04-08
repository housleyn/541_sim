import os 
import sys
import pytest 
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from mesh import Mesh
from boundaries import Boundary
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
    mesh = Mesh(domain, nx=20, ny=2)
    return mesh

@pytest.fixture
def boundary(mesh):
    bc = Boundary(mesh)
    return bc

def test_nodes_on_left_boundary_set_correctly(boundary):
    # Make sure mesh is constructed
    boundary.mesh.construct_mesh()

    # Now apply the boundary conditions
    boundary.set_left_boundary(var='p', value=10.0)
    boundary.apply()

    # Check that left boundary nodes have b = 10.0
    for j, i in boundary.mesh.left_node_indices:
        node = boundary.mesh.nodes[j, i]
        assert node.b == 10.0
        assert node.aP == 1.0


def test_nodes_on_right_boundary_set_correctly(boundary):
    
    boundary.set_right_boundary(var='p', value=20.0)
    boundary.apply()
    mesh = boundary.mesh
    for j in range(mesh.nodes.shape[0]):
        node = mesh.nodes[j, -1]
        if node is not None:
            assert node.aP == 1
            assert node.aE == 0
            assert node.aW == 0
            assert node.b == 20.0

def test_u_faces_on_left_boundary_set_correctly(boundary):
    
    boundary.set_left_boundary(var='u', value=1.0)
    boundary.apply()
    mesh = boundary.mesh
    for j in range(mesh.u_faces.shape[0]):
        if mesh.u_faces[j, 0] is not None:
            assert mesh.u_faces[j, 0].u == 1.0
        

def test_u_faces_on_right_boundary_set_correctly(boundary):
    
    boundary.set_right_boundary(var='u', value=0.5)
    boundary.apply()
    mesh = boundary.mesh
    for j in range(mesh.u_faces.shape[0]):
        if mesh.u_faces[j, -1] is not None:
            assert mesh.u_faces[j, -1].u == 0.5
        

def test_v_faces_on_upper_boundary_set_correctly(boundary):
    
    boundary.set_upper_boundary(var='v', value=0.0)
    boundary.apply()
    mesh = boundary.mesh
    for i in range(mesh.v_faces.shape[1]):
        if mesh.v_faces[-1, i] is not None:  
            assert mesh.v_faces[-1, i].v == 0.0

def test_v_faces_on_lower_boundary_set_correctly(boundary):
    
    boundary.set_lower_boundary(var='v', value=0.0)
    boundary.apply()
    mesh = boundary.mesh
    for i in range(mesh.v_faces.shape[1]):
        if mesh.v_faces[0, i] is not None:
            assert mesh.v_faces[0, i].v == 0.0



def test_1D_boundary_only_sets_left_and_right_nodes():
    d = Domain()
    d.define_left_boundary(lambda y, x=None: 0)
    d.define_right_boundary(lambda y, x=None: 1)
    d.define_lower_boundary(lambda x, y=None: 0)  # dummy flat bottom
    d.define_upper_boundary(lambda x, y=None: 0)  # dummy flat top
    mesh = Mesh(d, nx=20)
    mesh.construct_mesh()

    boundary = Boundary(mesh)
    boundary.set_left_boundary(var='p', value=5.0)
    boundary.set_right_boundary(var='p', value=15.0)
    boundary.apply()

    # Should only affect 1D nodes, no v_faces, u_faces involved
    for i in range(mesh.shape[0]):
        node = mesh.nodes[i]
        if node is not None:
            if i == 0:
                assert node.aP == 1
                assert node.b == 5.0
            elif i == mesh.shape[0] - 1:
                assert node.aP == 1
                assert node.b == 15.0
            else:
                # Interior nodes should remain untouched
                assert node.aP != 1 or node.b not in [5.0, 15.0]


def test_setting_only_pressure_does_not_set_velocity(boundary):
    boundary.set_left_boundary(var='p', value=10.0)
    boundary.apply()
    mesh = boundary.mesh

    for j in range(mesh.u_faces.shape[0]):
        if mesh.u_faces[j, 0] is not None:
            # Ensure velocity is still 0.0 or default
            assert mesh.u_faces[j, 0].u == 0.0
    for j in range(mesh.u_faces.shape[0]):
        if mesh.u_faces[j, -1] is not None:
            assert mesh.u_faces[j, -1].u == 0.0


def test_setting_only_velocity_does_not_set_pressure(boundary):
    boundary.set_right_boundary(var='u', value=1.0)
    boundary.apply()
    mesh = boundary.mesh

    for j in range(mesh.nodes.shape[0]):
        node = mesh.nodes[j, -1]
        if node is not None:
            # Ensure pressure coefficients are still unset
            assert node.aP != 1 or node.b != 1.0
