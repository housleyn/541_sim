from domain import Domain 
from mesh import Mesh
from material import Material
from boundaries import Boundary
from simple import SIMPLE


#construct domain
domain = Domain()
domain.define_lower_boundary(lambda x, y: y - .1 * x + 0.25)
domain.define_upper_boundary(lambda x, y: y + .1 * x - 0.25)
domain.define_left_boundary(lambda y, x=None:0)
domain.define_right_boundary(lambda y, x=None:2)


#construct mesh
mesh = Mesh(domain, nx=20)
mesh.construct_mesh()
# mesh.plot_mesh()


#define material properties
fluid = Material()
fluid.rho = 1.0
fluid.mu = 0

#define boundary conditions
boundary = Boundary(mesh)
boundary.set_left_boundary('p', 10)
boundary.set_right_boundary('p', 0)
boundary.apply()


#run simple class
simple = SIMPLE(mesh, boundary, fluid)
simple.run()
