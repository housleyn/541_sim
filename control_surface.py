class ControlSurface:
    def __init__(self):
        self.u = None
        self.v = None
        self.orientation = None  # Will be set when u or v is assigned
        self.area = None
        self.position = None
        self.b = None

    @property
    def u(self):
        return self._u

    @u.setter
    def u(self, value):
        self._u = value
        self.orientation = 'vertical'  # u is aligned with x-direction → vertical face

    @property
    def v(self):
        return self._v

    @v.setter
    def v(self, value):
        self._v = value
        self.orientation = 'horizontal'  # v is aligned with y-direction → horizontal face
