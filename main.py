from domain import Domain 
from mesh import Mesh



domain = Domain()
domain.define_lower_boundary(lambda x, y: y - .1 * x + 0.25)
domain.define_upper_boundary(lambda x, y: y + .1 * x - 0.25)
domain.define_left_boundary(lambda y, x=None:0)
domain.define_right_boundary(lambda y, x=None:2)
# domain.plot_domain(x_range=(0, 2), y_range=(-0.25, 0.25))


mesh = Mesh(domain, dx=0.1)
mesh.construct_nodes()
print("Mesh node positions:")
if len(mesh.shape) == 1:  # 1D mesh
    for i, node in enumerate(mesh.nodes):
        if node is not None:
            print(f"Node[{i}] → {node.position}")
else:  # 2D mesh
    for j in range(mesh.shape[0]):     # y-direction (rows)
        for i in range(mesh.shape[1]): # x-direction (columns)
            node = mesh.nodes[j, i]
            if node is not None:
                print(f"Node[{j}, {i}] → {node.position}")
mesh.plot_mesh()