from node import Node
from control_surface import ControlSurface
import numpy as np 
import matplotlib.pyplot as plt
import scipy
from scipy.sparse import lil_matrix, csr_matrix

class Mesh():
    def __init__(self, domain, nx, ny=1):
        self.domain = domain
        self.nx = nx
        self.ny = ny
        
        
        
        self.alpha_u = None
        self.alpha_v = None
        
        x_bounds, y_bounds = self.domain.get_bounds()
        x_min, x_max = x_bounds
        y_min, y_max = y_bounds

        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max

        self.dx = np.full(nx, (x_max - x_min) / nx)
        self.x_centers = np.linspace(x_min + 0.5 * self.dx[0], x_max - 0.5 * self.dx[0], nx)

        if ny == 1:
            self.dy = np.array([
                self.domain.upper_boundary_func(x, 0) - self.domain.lower_boundary_func(x, 0)
                for x in self.x_centers
            ])
            self.y_centers = np.zeros(nx)
            self.shape = (nx,)
        else:
            self.dy = np.full(ny, (y_max - y_min) / ny)
            self.y_centers = np.linspace(y_min + 0.5 * self.dy[0], y_max - 0.5 * self.dy[0], ny)
            self.shape = (ny, nx)



        self.nodes = np.empty(self.shape, dtype=object)
        if ny > 1:
            self.u_faces = np.empty((ny, nx - 1), dtype=object)
            self.v_faces = np.empty((ny - 1, nx), dtype=object)
        else:
            self.u_faces = np.empty((1, nx - 1), dtype=object)
            self.v_faces = None

        self.build_valid_index_lists()

    def construct_nodes(self):
        if self.ny == 1:
            for i in range(self.nx):
                x = self.x_min + sum(self.dx[:i])
                y = 0.0
                if self.domain.is_inside(x, y):
                    node = Node()
                    node.position = (x, y)
                    self.nodes[i] = node
                else:
                    self.nodes[i] = None
        else:
            for j in range(self.ny):
                for i in range(self.nx):
                    x = self.x_min + sum(self.dx[:i])
                    y = self.y_min + sum(self.dy[:j])
                    if self.domain.is_inside(x, y):
                        node = Node()
                        node.position = (x, y)
                        self.nodes[j, i] = node
                    else:
                        self.nodes[j, i] = None

    def construct_control_surfaces(self):
        is_2d = len(self.shape) == 2
        nx = self.shape[1] if is_2d else self.shape[0]
        ny = self.shape[0] if is_2d else 1

        # For u-surfaces (control surfaces between nodes in x-direction)
        for j in range(ny):
            for i in range(nx - 1):  # We have nx - 1 u-faces for each row
                n1 = self.nodes[j, i] if is_2d else self.nodes[i]
                n2 = self.nodes[j, i + 1] if is_2d else self.nodes[i + 1]

                # Ensure that both nodes are not None
                if n1 is None or n2 is None:
                    self.u_faces[j, i] = None if is_2d else None
                    continue

                # Create a control surface at the average of the node positions
                surf = ControlSurface()
                x = 0.5 * (n1.position[0] + n2.position[0])
                y = 0.5 * (n1.position[1] + n2.position[1])
                surf.position = (x, y)
                surf.u = 0.0  # This triggers the orientation to be set as 'vertical'
                self.u_faces[j, i] = surf if is_2d else None

        # For v-surfaces (control surfaces between nodes in y-direction)
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
                    surf.v = 0.0  # This triggers the orientation to be set as 'horizontal'
                    self.v_faces[j, i] = surf





    def build_A_matrix(self, var_type):
        is_2d = len(self.shape) == 2
        if var_type == 'u':
            face_array = self.u_faces
            indices = list(self.u_face_indices)
        elif var_type == 'v':
            if not is_2d:
                return None, None
            face_array = self.v_faces
            indices = list(self.v_face_indices)
        elif var_type == 'p':
            face_array = self.nodes
            indices = self.interior_node_indices
        else:
            raise ValueError(f"Unknown variable type: {var_type}")

        N = len(indices)
        A = lil_matrix((N, N))
        b = np.zeros(N)
        index_map = {idx: k for k, idx in enumerate(indices)}

        for k, idx in enumerate(indices):
            if is_2d:
                obj = face_array[idx[0], idx[1]]
            else:
                obj = face_array[idx]
            if obj is None:
                continue

            if hasattr(obj, "bc_type") and obj.bc_type == "Dirichlet":
                A[k, k] = 1.0
                if var_type == 'u':
                    b[k] = obj.u if obj.u is not None else 0.0
                elif var_type == 'v':
                    b[k] = obj.v if obj.v is not None else 0.0
                elif var_type == 'p':
                    b[k] = obj.p if obj.p is not None else 0.0
                continue

            relax = self.alpha_u if var_type == 'u' else self.alpha_v if var_type == 'v' else 1.0
            A[k, k] = (obj.aP / relax) if obj.aP is not None else 1e-12
            b[k] = obj.b if obj.b is not None else 0.0

            if is_2d:
                neighbors = {
                    'aE': (idx[0], idx[1] + 1),
                    'aW': (idx[0], idx[1] - 1),
                    'aN': (idx[0] + 1, idx[1]),
                    'aS': (idx[0] - 1, idx[1])
                }
            else:
                neighbors = {
                    'aE': idx + 1,
                    'aW': idx - 1
                }

            for coeff_name, neighbor_idx in neighbors.items():
                coeff = getattr(obj, coeff_name, None)
                if coeff is not None and neighbor_idx in index_map:
                    A[k, index_map[neighbor_idx]] = -coeff

        return csr_matrix(A), b



    def build_valid_index_lists(self):
        self.interior_node_indices = []
        self.u_face_indices = []
        self.v_face_indices = []

        self.left_node_indices = []
        self.right_node_indices = []
        self.top_node_indices = []
        self.bottom_node_indices = []

        if self.ny == 1:
            # 1D case
            for i in range(self.nx):
                if self.nodes[i] is not None:
                    self.interior_node_indices.append(i)
                    if i == 0:
                        self.left_node_indices.append(i)
                    elif i == self.nx - 1:
                        self.right_node_indices.append(i)

            for i in range(self.nx - 1):
                if self.u_faces[0, i] is not None:
                    self.u_face_indices.append((0, i))

        else:
            # 2D case
            for j in range(self.ny):
                for i in range(self.nx):
                    if self.nodes[j, i] is not None:
                        self.interior_node_indices.append((j, i))

                        if i == 0:
                            self.left_node_indices.append((j, i))
                        elif i == self.nx - 1:
                            self.right_node_indices.append((j, i))

                        if j == 0:
                            self.bottom_node_indices.append((j, i))
                        elif j == self.ny - 1:
                            self.top_node_indices.append((j, i))

            for j in range(self.ny):
                for i in range(self.nx - 1):
                    if self.u_faces[j, i] is not None:
                        self.u_face_indices.append((j, i))

            for j in range(self.ny - 1):
                for i in range(self.nx):
                    if self.v_faces[j, i] is not None:
                        self.v_face_indices.append((j, i))


    
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

    
