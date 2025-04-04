from node import Node
from control_surface import ControlSurface
import numpy as np 
import matplotlib.pyplot as plt

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

    
