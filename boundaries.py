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
        self.apply_pressure_boundary('left', self.left_boundary)
        self.apply_pressure_boundary('right', self.right_boundary)
        self.apply_pressure_boundary('top', self.upper_boundary)
        self.apply_pressure_boundary('bottom', self.lower_boundary)

        self.apply_velocity_boundary('left', 'u', self.left_boundary)
        self.apply_velocity_boundary('right', 'u', self.right_boundary)
        self.apply_velocity_boundary('top', 'v', self.upper_boundary)
        self.apply_velocity_boundary('bottom', 'v', self.lower_boundary)

    def apply_pressure_boundary(self, side, bc):
        if bc and bc['var'] == 'p':
            indices = self.mesh.get_boundary_node_indices(side)

            if len(self.mesh.shape) == 2:
                # 2D: unpack j, i
                for j, i in indices:
                    node = self.mesh.nodes[j, i]
                    node.aP = 1
                    if side in ['left', 'right']:
                        node.aW = 0
                        node.aE = 0
                    else:
                        node.aN = 0
                        node.aS = 0
                    node.b = bc['value']
            else:
                # 1D: just i
                for i in indices:
                    node = self.mesh.nodes[i]
                    node.aP = 1
                    node.aW = 0
                    node.aE = 0
                    node.b = bc['value']


    def apply_velocity_boundary(self, side, var, bc):
        if bc and bc['var'] == var:
            if var == 'u':
                face_array = self.mesh.u_faces
                idx = 0 if side == 'left' else -1
                for j in range(face_array.shape[0]):
                    face = face_array[j, idx]
                    if face is not None:
                        face.u = bc['value']
            elif var == 'v':
                face_array = self.mesh.v_faces
                idx = 0 if side == 'bottom' else -1
                for i in range(face_array.shape[1]):
                    face = face_array[idx, i]
                    if face is not None:
                        face.v = bc['value']
