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
mesh.plot_mesh()