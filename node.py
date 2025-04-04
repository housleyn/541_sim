class Node():
    def __init__(self):
        self.aW = None
        self.aE = None
        self.aP = None 
        self.aN = None 
        self.aS = None 
        self.b = None 
        self.P = None
        self.position = None 

    def define_coefficients(self, De, Fe, Dw, Fw, Dn, Fn, Ds, Fs, Sc, Sp, dx, dy):
        self.aE = De + max(-Fe,0)
        self.aW = Dw + max(Fw,0)        
        self.aN = Dn + max(-Fn,0)
        self.aS = Ds + max(Fs,0)
        self.b = Sc * dx* dy
        self.aP = self.aW + self.aE + self.aN + self.aS - Sp*dx*dy 
        
        