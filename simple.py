import numpy as np

class SIMPLE():
    def __init__(self, mesh, material):
        self.rho = material.rho
        self.mu = material.mu
        self.mesh = mesh
        self.domain = mesh.domain
        is_2d = len(mesh.shape) == 2
        ny, nx = mesh.shape if is_2d else (1, mesh.shape[0])


        self.u_star_old = np.zeros((ny, nx-1))
        self.u_star_new = np.zeros_like(self.u_star_old)
        self.v_star_old = np.zeros((ny-1, nx)) if is_2d else None
        self.v_star_new = np.zeros_like(self.v_star_old) if is_2d else None
        
        self.p_star_old = np.zeros((ny, nx))
        self.p_star_new = np.zeros_like(self.p_star_old)
        self.p_prime = np.zeros_like(self.p_star_old)

        self.dx = mesh.dx
        self.dy = mesh.dy if is_2d else 1.0
        self.alph_u = .7
        self.alpha_v = .3 

    def calculate_diffusion_terms(self):
        De = Dw = self.mu / self.dx 
        Dn = Ds = self.mu / self.dy if len(self.mesh.shape) == 2 else 0.0
        return De, Dw, Dn, Ds 
    
    def generate_initial_guesses(self):
        self.u_star_old = 
        if self.v_star_old is not None:
            self.v_star_old =
        self.p_star_old =
    
    def solve_discretized_x_momentum(self):
        dx, dy , rho, mu = self.dx, self.dy, self.rho, self.mu
        De, Dw, Dn, Ds = self.calculate_diffusion_terms()
        for j, i in self.mesh.u_face_indices:
            Fe = (self.rho/2) * (self.u_star_old[j,i+1] + self.u_star_old[j,i])
            Fw = (self.rho/2) * (self.u_star_old[j,i] + self.u_star_old[j,i-1])
            Fn = (self.rho/2) * (self.v_star_old[j+1,i] + self.v_star_old[j+1,i-1]) if self.v_star_old is not None else 0.0
            Fs = (self.rho/2) * (self.v_star_old[j,i] + self.v_star_old[j,i-1]) if self.v_star_old is not None else 0.0
            #I'm stuck, I need to read all the files.
            

        
    
    def solve_discretized_y_momentum(self):
        self.v_star_new =

    def solve_pressure_correction(self):
        self.p_prime =

    def calculate_new_pressure_field(self):
        self.p_new =

    def calcualte_new_velocity_field(self):
        self.u_star_new =
        if self.v_star_old is not None:
            self.v_star_new =

    def set_old_equal_to_new(self):
        self.u_star_old = self.u_star_new
        if self.v_star_old is not None:
            self.v_star_old = self.v_star_new
        self.p_star_old = self.p_new
    
    


    def run(self): #iterate until convergence
        is_2d = len(self.mesh.shape) == 2
        self.generate_initial_guesses()
        converged = False

        while not converged:
            self.solve_discretized_x_momentum()
            if is_2d:
                self.solve_discretized_y_momentum()
            self.solve_pressure_correction()
            self.calculate_new_pressure_field()
            self.calcualte_new_velocity_field()
            self.set_old_equal_to_new()