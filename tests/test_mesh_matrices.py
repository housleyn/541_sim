import os 
import sys 
import pytest
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from mesh import Mesh
from scipy.sparse import csr_matrix
from mesh import Mesh
from domain import Domain


@pytest.fixture
def domain():
    d = Domain()
    d.define_lower_boundary(lambda x, y: y - .1 * x + 0.25)
    d.define_upper_boundary(lambda x, y: y + .1 * x - 0.25)
    d.define_left_boundary(lambda y, x=None: 0)
    d.define_right_boundary(lambda y, x=None: 2)
    return d

@pytest.fixture
def mesh(domain):
    m = Mesh(domain, nx=20, ny=2)
    return m


def test_matrix_build_var_type_switch(mesh):
    for var in ['u', 'v', 'p']:
        A, b = mesh.build_A_matrix(var)
        assert isinstance(A, csr_matrix), f"A should be a sparse matrix for {var}"
        assert b.ndim == 1, f"b should be a 1D vector for {var}"


def test_matrix_shapes_match(mesh):
    for var in ['u', 'v', 'p']:
        A, b = mesh.build_A_matrix(var)
        assert A.shape[0] == A.shape[1], "A must be square"
        assert A.shape[0] == b.shape[0], "A and b must match in length"


def test_boundary_conditions_applied(mesh):
    A, b = mesh.build_A_matrix('u')
    index_map = {idx: k for k, idx in enumerate(mesh.u_face_indices)}
    for idx in mesh.u_face_indices:
        face = mesh.u_faces[idx]
        if hasattr(face, "bc_type") and face.bc_type == "Dirichlet":
            row = index_map[idx]
            assert A[row, row] == pytest.approx(1.0)
            assert b[row] == pytest.approx(face.u or 0.0)


def test_matrix_coefficients_inserted_correctly(mesh):
    A, b = mesh.build_A_matrix('u')
    index_map = {idx: k for k, idx in enumerate(mesh.u_face_indices)}
    for idx, k in index_map.items():
        face = mesh.u_faces[idx]
        if hasattr(face, "bc_type") and face.bc_type == "Dirichlet":
            continue  # Already checked
        assert A[k, k] == pytest.approx(face.aP or 1e-12)
        assert b[k] == pytest.approx(face.b or 0.0)
        neighbors = {
            'aE': (idx[0], idx[1] + 1),
            'aW': (idx[0], idx[1] - 1),
            'aN': (idx[0] + 1, idx[1]),
            'aS': (idx[0] - 1, idx[1])
        }
        for coeff_name, neighbor_idx in neighbors.items():
            if neighbor_idx in index_map:
                coeff = getattr(face, coeff_name, None)
                if coeff is not None:
                    col = index_map[neighbor_idx]
                    assert A[k, col] == pytest.approx(-coeff)