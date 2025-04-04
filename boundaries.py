class Boundary():
    def __init__(self, mesh):
        self.mesh = mesh
        self.left_boundary = None
        self.right_boundary = None
        self.upper_boundary = None
        self.lower_boundary = None

    def set_left_boundary(self, var, value, kind='Dirichlet'):
        self.left_boundary = {'var': var, 'value': value, 'kind': kind}

    def set_right_boundary(self, var, value, kind='Dirichlet'):
        self.right_boundary = {'var': var, 'value': value, 'kind': kind}

    def set_upper_boundary(self, var, value, kind='Dirichlet'):
        self.upper_boundary = {'var': var, 'value': value, 'kind': kind}

    def set_lower_boundary(self, var, value, kind='Dirichlet'):
        self.lower_boundary = {'var': var, 'value': value, 'kind': kind}

    def apply(self):
        
        nodes = self.mesh.nodes
        is_2d = len(nodes.shape) == 2
        u_faces = self.mesh.u_faces
        v_faces = self.mesh.v_faces

        # LEFT boundary
        if self.left_boundary:
            var = self.left_boundary['var']
            value = self.left_boundary['value']
            if var == 'p':  # pressure stored at nodes
                for j in range(nodes.shape[0]):
                    node = nodes[j,0]
                    if node is not None:
                        node.aP = 1
                        node.aW = 0
                        node.aE = 0
                        node.b = value
            elif var == 'u':  # u velocity stored at control surface
                for j in range(u_faces.shape[0]):
                    if u_faces[j,0] is not None:
                        u_faces[j,0].u = value
                

        # RIGHT boundary
        if self.right_boundary:
            var = self.right_boundary['var']
            value = self.right_boundary['value']
            if var == 'p':
                for j in range(nodes.shape[0]):
                    node = nodes[j,-1]
                    if node is not None:
                        node.aP = 1
                        node.aW = 0
                        node.aE = 0
                        node.b = value
            elif var == 'u':
                for j in range(u_faces.shape[0]):
                    if u_faces[j,-1] is not None:
                        u_faces[j,-1].u = value

        # LOWER boundary (j=0 for all i)
        if self.lower_boundary:
            var, value = self.lower_boundary['var'], self.lower_boundary['value']
            if var == 'p' and is_2d:
                for i in range(nodes.shape[1]):
                    node = nodes[0, i]
                    if node is not None:
                        node.aP = 1; node.aN = 0; node.aS = 0; node.b = value
            elif var == 'v':
                for i in range(v_faces.shape[1]):
                    if v_faces[0, i] is not None:
                        v_faces[0, i].v = value

        # UPPER boundary (j=ny-1 for all i)
        if self.upper_boundary:
            var, value = self.upper_boundary['var'], self.upper_boundary['value']
            if var == 'p' and is_2d:
                for i in range(nodes.shape[1]):
                    node = nodes[-1, i]
                    if node is not None:
                        node.aP = 1; node.aN = 0; node.aS = 0; node.b = value
            elif var == 'v':
                for i in range(v_faces.shape[1]):
                    if v_faces[-1, i] is not None:
                        v_faces[-1, i].v = value
