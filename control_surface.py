class ControlSurface:
    def __init__(self):
        self.u = None
        self.v = None
        self.orientation = None  # Will be set when u or v is assigned
        self.area = None
        self.position = None
        self.b = None
        self.aE = None 
        self.aW = None
        self.aN = None
        self.aS = None
        self.aP = None
        self.u_old = None
        self.v_old = None

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


  
    def calculate_x_coefficients(self, dx, dy, De, Fe, Dw, Fw, Dn, Fn, Ds, Fs,pe, pw, alphau,  is_2d):
        self.aE = De * dy + max(-Fe, 0)*dy 
        self.aW = Dw*dy + max(Fw,0)*dy 
        self.aN = Dn*dx + max(-Fn,0)*dx if is_2d else None
        self.aS = Ds*dx + max(Fs,0)*dx if is_2d else None
        self.aP = self.aE + self.aW + self.aN + self.aS + (Fe-Fw)*dy + (Fn-Fs)*dx 
        self.b = (pe-pw)*dy + ((1-alphau)*self.aP/alphau)*self.u_old

    def calculate_y_coefficients(self, dx, dy, De, Fe, Dw, Fw, Dn, Fn, Ds, Fs, ps, pn, alphav, is_2d):
        self.aE = De * dy + max(-Fe, 0)*dy if is_2d else None
        self.aW = Dw*dy + max(Fw,0)*dy if is_2d else None
        self.aN = Dn*dx + max(-Fn,0)*dx 
        self.aS = Ds*dx + max(Fs,0)*dx
        self.aP = self.aE + self.aW + self.aN + self.aS + (Fe-Fw)*dy + (Fn-Fs)*dx 
        self.b = (ps-pn)*dx + ((1-alphav)*self.aP/alphav)*self.v_old
