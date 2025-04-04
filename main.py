from domain import Domain 
from mesh import Mesh
from material import Material


#construct domain
domain = Domain()
domain.define_lower_boundary(lambda x, y: y - .1 * x + 0.25)
domain.define_upper_boundary(lambda x, y: y + .1 * x - 0.25)
domain.define_left_boundary(lambda y, x=None:0)
domain.define_right_boundary(lambda y, x=None:2)


#construct mesh
mesh = Mesh(domain, dx=0.1)
mesh.construct_nodes()
mesh.construct_control_surfaces()
# mesh.plot_mesh()


#define material properties
fluid = Material()
fluid.rho = 1.0
fluid.mu = 0

#define boundary conditions



#run simple class