class Node():
    def __init__(self):
        self.aW = None
        self.aE = None
        self.aP = None 
        self.aN = None 
        self.aS = None 
        self.b = None 
        self.p = None
        self.position = None 

    def define_pressure_correction_coefficients(self, alphau, alphav,rho, dx, dy,
                                                aP_u_E, aP_u_W,
                                                aP_v_N, aP_v_S,
                                                u_star_W, u_star_E,
                                                v_star_S, v_star_N):
        self.aE = rho * dy**2 * alphau / aP_u_E if aP_u_E else 0.0
        self.aW = rho * dy**2 * alphau / aP_u_W if aP_u_W else 0.0
        self.aN = rho * dx**2 * alphav / aP_v_N if aP_v_N else 0.0
        self.aS = rho * dx**2 * alphav / aP_v_S if aP_v_S else 0.0

        self.aP = self.aE + self.aW + self.aN + self.aS

        self.b = rho * ((u_star_W - u_star_E) * dy + (v_star_S - v_star_N) * dx)

        
            
            