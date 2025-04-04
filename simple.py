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


    def generate_initial_guesses(self):
        self.u_star_old = 
        if self.v_star_old is not None:
            self.v_star_old =
        self.p_star_old =
    
    def solve_discretized_x_momentum(self):
        is_2d = len(self.mesh.shape) == 2
        ny, nx = self.mesh.shape if is_2d else (1, self.mesh.shape[0])
        for j in range(ny):
            for i in range(nx-1):
                Fs = #part 3 notes

        self.u_star_new = 
    
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