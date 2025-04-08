import os 
import sys
import pytest
import numpy as np

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
    mesh = Mesh(domain, nx=20, ny=3)
    mesh.construct_mesh()
    return mesh
@pytest.fixture
def boundary(mesh):
    bc = Boundary(mesh)
    bc.set_left_boundary(var= 'p', value= 10)
    bc.set_right_boundary(var= 'p', value= 0)
    bc.set_upper_boundary(var= 'v', value= 0, kind='Dirichlet')
    bc.set_lower_boundary(var= 'u', value= 0, kind='Dirichlet')
    bc.apply()
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

def test_xmin_not_equal_xmax(simple):
    xmin = simple.mesh.domain.x_min
    xmax = simple.mesh.domain.x_max
    print(f"xmin = {xmin}, xmax = {xmax}")
    assert xmin != xmax, "Domain x_min and x_max are equal — interpolation will fail."

def test_interior_node_indices_not_empty(simple):
    simple.mesh.construct_nodes()
    simple.mesh.construct_control_surfaces()
    simple.mesh.build_valid_index_lists()

    assert len(simple.mesh.interior_node_indices) > 0, "interior_node_indices is empty, but should contain valid node indices."

def test_all_nodes_have_information_set(simple):
    # You must trigger the logic that sets node coefficients
    simple.generate_initial_guesses()
    simple.solve_discretized_x_momentum()
    if simple.v_star_old is not None:
        simple.solve_discretized_y_momentum()
    simple.solve_pressure_correction()

    for idx in simple.mesh.interior_node_indices:
        if len(simple.mesh.shape) == 2:
            j, i = idx
            node = simple.mesh.nodes[j, i]
        else:
            i = idx
            j = 0
            node = simple.mesh.nodes[i]

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

def test_simple_initialization(mesh,boundary, material):
    s = SIMPLE(mesh, boundary,material)
    assert np.all(s.u_star_old == 0)
    assert s.v_star_old is None or np.all(s.v_star_old == 0)
    assert s.p_star_old.shape == (1, mesh.shape[0]) if len(mesh.shape) == 1 else mesh.shape

def test_calculate_diffusion_terms_1d(simple):
    simple.mesh.shape = (10,)  # force 1D
    De, Dw, Dn, Ds = simple.calculate_diffusion_terms()

    # Ensure that all values in De and Dw match simple.mu / simple.dx
    expected_value = simple.mu / simple.dx
    assert np.allclose(De, expected_value), f"De values: {De} do not match expected value: {expected_value}"
    assert np.allclose(Dw, expected_value), f"Dw values: {Dw} do not match expected value: {expected_value}"


def test_calculate_diffusion_terms_2d(simple):
    De, Dw, Dn, Ds = simple.calculate_diffusion_terms()

    # Ensure that all values in De and Dw match simple.mu / simple.dx
    expected_value = simple.mu / simple.dx
    assert np.allclose(De, expected_value), f"De values: {De} do not match expected value: {expected_value}"
    assert np.allclose(Dw, expected_value), f"Dw values: {Dw} do not match expected value: {expected_value}"


def test_generate_initial_guesses_pressure(simple):
    # Force pressure boundary conditions

    simple.generate_initial_guesses()
    assert not np.all(simple.p_star_old == 0)

def test_generate_initial_guesses_velocity_u(simple):
    simple.boundary.left_boundary = {'var': 'u', 'value': 3.0}
    simple.generate_initial_guesses()
    assert np.all(simple.u_star_old == 3.0)

def test_generate_initial_guesses_velocity_v(simple):
    simple.boundary.lower_boundary = {'var': 'v', 'value': 2.0}
    simple.generate_initial_guesses()
    if simple.v_star_old is not None:
        assert np.all(simple.v_star_old == 2.0)

def test_solve_discretized_x_momentum(simple):
    simple.generate_initial_guesses()
    simple.solve_discretized_x_momentum()
    assert not np.all(simple.u_star_new == simple.u_star_old)

def test_solve_discretized_y_momentum(simple):
    if simple.v_star_old is not None:
        simple.generate_initial_guesses()
        simple.solve_discretized_x_momentum()
        simple.solve_discretized_y_momentum()
        assert not np.all(simple.v_star_new == simple.v_star_old)

def test_solve_pressure_correction(simple):
    simple.generate_initial_guesses()
    simple.solve_discretized_x_momentum()
    if simple.v_star_old is not None:
        simple.solve_discretized_y_momentum()
    simple.solve_pressure_correction()
    assert not np.all(simple.p_prime == 0)

def test_calculate_new_pressure_field(simple):
    simple.p_star_old.fill(1.0)
    simple.p_prime.fill(0.5)
    simple.alph_p = 0.7
    simple.calculate_new_pressure_field()
    assert np.any(simple.p_star_new != 1.0)

def test_calculate_new_velocity_field(simple):
    simple.u_star_new.fill(1.0)
    simple.p_prime.fill(0.5)
    simple.calculate_new_velocity_field()
    assert np.any(simple.u_star_new != 1.0)

def test_set_old_equal_to_new(simple):
    simple.u_star_new.fill(9.0)
    if simple.v_star_new is not None:
        simple.v_star_new.fill(7.0)
    simple.p_star_new.fill(5.0)
    simple.set_old_equal_to_new()
    assert np.all(simple.u_star_old == 9.0)
    if simple.v_star_old is not None:
        assert np.all(simple.v_star_old == 7.0)
    assert np.all(simple.p_star_old == 5.0)

def test_run_converges(simple):
    simple.tolerance = 1e5  # force early convergence
    simple.run()
    assert simple.iterations == 1
    assert np.linalg.norm(simple.p_prime) < simple.tolerance
