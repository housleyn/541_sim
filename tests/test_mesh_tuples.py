import os 
import sys 
import pytest 
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from mesh import Mesh
from domain import Domain
import pytest
import numpy as np

# Assuming you have access to Mesh and Domain classes
@pytest.fixture
def setup_2d_mesh():
    domain = Domain()
    domain.define_left_boundary(lambda y, x=None: 0)
    domain.define_right_boundary(lambda y, x=None: 2)
    domain.define_lower_boundary(lambda x, y: y - 0.1 * x + 0.25)
    domain.define_upper_boundary(lambda x, y: y + 0.1 * x - 0.25)
    
    mesh = Mesh(domain, nx=20,ny=3)
    mesh.construct_nodes()
    mesh.build_valid_index_lists()
    
    return mesh

@pytest.fixture
def setup_1d_mesh():
    domain = Domain()
    domain.define_left_boundary(lambda y, x=None: 0)
    domain.define_right_boundary(lambda y, x=None: 2)
    domain.define_lower_boundary(lambda x, y=0: 0)  # Lower boundary (no variation in y)
    domain.define_upper_boundary(lambda x, y=0: 0)  # Upper boundary (no variation in y)
    
    mesh = Mesh(domain, nx=20)
    mesh.construct_nodes()
    mesh.build_valid_index_lists()
    
    return mesh

# Test for 2D mesh
def test_interior_node_indices_2d(setup_2d_mesh):
    mesh = setup_2d_mesh
    # Check if interior_node_indices contains tuples for 2D
    for idx in mesh.interior_node_indices:
        assert isinstance(idx, tuple), f"Expected tuple, got {type(idx)}"
        assert len(idx) == 2, f"Expected tuple of length 2, got {len(idx)}"

# Test for 1D mesh
def test_interior_node_indices_1d(setup_1d_mesh):
    mesh = setup_1d_mesh
    # Check if interior_node_indices contains only single indices for 1D
    for idx in mesh.interior_node_indices:
        assert isinstance(idx, int), f"Expected int, got {type(idx)}"
