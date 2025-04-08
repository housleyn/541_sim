from node import Node
from control_surface import ControlSurface
import numpy as np 
import matplotlib.pyplot as plt
import scipy
from scipy.sparse import lil_matrix, csr_matrix

class Mesh():
    def __init__(self, domain, dx=None, dy=None):
        self.domain = domain
        self.dx = dx 
        self.dy = dy
        self.shape = np.zeros
        
        x_bounds, y_bounds = self.domain.get_bounds()
        x_min, x_max = x_bounds
        y_min, y_max = y_bounds

        nx = int((x_max - x_min) / dx) + 1 if dx else 1
        ny = int((y_max - y_min) / dy) + 1 if dy else 1

        self.shape = (ny, nx) if dy else (nx,)  # 1D or 2D        
        self.nodes = np.empty(self.shape, dtype=Node)

        self.u_faces = np.empty((ny, nx - 1), dtype=ControlSurface)
        self.v_faces = np.empty((ny - 1, nx), dtype=ControlSurface)

        self.build_valid_index_lists()

    def construct_nodes(self):
        x_bounds, y_bounds = self.domain.get_bounds()
        x_min, x_max = x_bounds
        y_min, y_max = y_bounds

        is_2d = len(self.shape) == 2
        nx = self.shape[1] if is_2d else self.shape[0]
        ny = self.shape[0] if is_2d else 1

        for j in range(ny):
            for i in range(nx):
                x = x_min + i * self.dx
                if is_2d:
                    y = (y_min + j * self.dy) 
                else:
                    y = 0.0
                if self.domain.is_inside(x, y):
                    node = Node()
                    node.position = (x, y)
                    if is_2d:
                        self.nodes[j, i] = node
                    else:
                        self.nodes[i] = node
                else:
                    if is_2d:
                        self.nodes[j, i] = None
                    else:
                        self.nodes[i] = None

    def construct_control_surfaces(self):
        is_2d = len(self.shape) == 2
        nx = self.shape[1] if is_2d else self.shape[0]
        ny = self.shape[0] if is_2d else 1


        x_bounds, y_bounds = self.domain.get_bounds()
        x_min, y_min = x_bounds[0], y_bounds[0]

        # u-surfaces (between nodes in x)
        for j in range(ny):
            for i in range(nx - 1):
                n1 = self.nodes[j, i] if is_2d else self.nodes[i]
                n2 = self.nodes[j, i + 1] if is_2d else self.nodes[i + 1]
                if n1 is None or n2 is None:
                    self.u_faces[j, i] = None if is_2d else None
                    continue
                surf = ControlSurface()
                x = 0.5 * (n1.position[0] + n2.position[0])
                y = 0.5 * (n1.position[1] + n2.position[1])
                surf.position = (x, y)
                surf.u = 0.0  # triggers orientation
                self.u_faces[j, i] = surf if is_2d else None

        # v-surfaces (between nodes in y)
        if is_2d:
            for j in range(ny - 1):
                for i in range(nx):
                    n1 = self.nodes[j, i]
                    n2 = self.nodes[j + 1, i]
                    if n1 is None or n2 is None:
                        self.v_faces[j, i] = None
                        continue
                    surf = ControlSurface()
                    x = 0.5 * (n1.position[0] + n2.position[0])
                    y = 0.5 * (n1.position[1] + n2.position[1])
                    surf.position = (x, y)
                    surf.v = 0.0  # triggers orientation
                    self.v_faces[j, i] = surf



    def build_A_matrix(self, var_type):
        if var_type == 'u':
            face_array = self.u_faces
            indices = list(self.u_face_indices)
        elif var_type == 'v':
            if len(self.shape) == 1:
                return None, None
            face_array = self.v_faces
            indices = list(self.v_face_indices)
        elif var_type == 'p':
            face_array = self.nodes
            indices = list(np.ndindex(self.nodes.shape))
        else:
            raise ValueError(f"Unknown variable type: {var_type}")

        N = len(indices)
        A = lil_matrix((N, N))
        b = np.zeros(N)
        index_map = {idx: k for k, idx in enumerate(indices)}

        for k, idx in enumerate(indices):
            obj = face_array[idx]
            if obj is None:
                continue

            # Handle Dirichlet BC directly
            if hasattr(obj, "bc_type") and obj.bc_type == "Dirichlet":
                A[k, k] = 1.0
                if var_type == 'u':
                    b[k] = obj.u if obj.u is not None else 0.0
                elif var_type == 'v':
                    b[k] = obj.v if obj.v is not None else 0.0
                elif var_type == 'p':
                    b[k] = obj.p if obj.p is not None else 0.0
                continue

            A[k, k] = obj.aP if obj.aP is not None else 1e-12
            b[k] = obj.b if obj.b is not None else 0.0

            neighbors = {
                'aE': (idx[0], idx[1] + 1),
                'aW': (idx[0], idx[1] - 1),
                'aN': (idx[0] + 1, idx[1]),
                'aS': (idx[0] - 1, idx[1])
            }

            for coeff_name, neighbor_idx in neighbors.items():
                coeff = getattr(obj, coeff_name, None)
                if coeff is not None and neighbor_idx in index_map:
                    A[k, index_map[neighbor_idx]] = -coeff

        return csr_matrix(A), b



    def build_valid_index_lists(self):
        self.interior_node_indices = []
        self.interior_nodes = []

        self.u_face_indices = []
        self.v_face_indices = []

        self.left_node_indices = []
        self.right_node_indices = []
        self.top_node_indices = []
        self.bottom_node_indices = []

        if hasattr(self, "nodes"):
            shape = self.nodes.shape

            if len(shape) == 2:
                ny, nx = shape
                for j in range(ny):
                    for i in range(nx):
                        node = self.nodes[j, i]
                        if node is not None:
                            self.interior_node_indices.append((j, i))
                            self.interior_nodes.append(node)

                            if i == 0:
                                self.left_node_indices.append((j, i))
                            elif i == nx - 1:
                                self.right_node_indices.append((j, i))
                            if j == 0:
                                self.bottom_node_indices.append((j, i))
                            elif j == ny - 1:
                                self.top_node_indices.append((j, i))

            elif len(shape) == 1:
                nx = shape[0]
                for i in range(nx):
                    node = self.nodes[i]
                    if node is not None:
                        self.interior_node_indices.append(i)
                        self.interior_nodes.append(node)
                        if i == 0:
                            self.left_node_indices.append(i)
                        elif i == nx - 1:
                            self.right_node_indices.append(i)

    
    def get_boundary_node_indices(self, side):
        if side == 'left':
            return self.left_node_indices
        elif side == 'right':
            return self.right_node_indices
        elif side == 'top':
            return self.top_node_indices
        elif side == 'bottom':
            return self.bottom_node_indices
        else:
            raise ValueError(f"Unknown boundary side: {side}")
    
    def construct_mesh(self):
        self.construct_nodes()
        self.construct_control_surfaces()
        self.build_valid_index_lists()

    def plot_mesh(self):
        x_vals = []
        y_vals = []

        if len(self.shape) == 1:
            for node in self.nodes:
                if node is not None:
                    x, y = node.position
                    x_vals.append(x)
                    y_vals.append(y)
        else:
            for row in self.nodes:
                for node in np.atleast_1d(row):
                    if node is not None:
                        x, y = node.position
                        x_vals.append(x)
                        y_vals.append(y)
        
        # Plot domain boundary
        x_bounds, y_bounds = self.domain.get_bounds()
        fig, ax = plt.subplots()
        self.domain.plot_domain(x_bounds, y_bounds)

        ax.scatter(x_vals, y_vals, s=10, c='red', label='Mesh Nodes')
        ax.legend()
        plt.show()

    
