import numpy as np
import matplotlib.pyplot as plt

class SIMPLE():
    def __init__(self, mesh, boundary, material):
        self.rho = material.rho
        self.mu = material.mu
        self.mesh = mesh
        self.domain = mesh.domain
        is_2d = len(mesh.shape) == 2
        ny, nx = mesh.shape if is_2d else (1, mesh.shape[0])
        self.iterations = 0
        self.tolerance = 1e-6
        self.boundary = boundary

        self.u_star_old = np.zeros((ny, nx - 1))
        self.u_star_new = np.zeros_like(self.u_star_old)
        self.v_star_old = np.zeros((ny - 1, nx)) if is_2d else None
        self.v_star_new = np.zeros_like(self.v_star_old) if is_2d else None

        self.p_star_old = np.zeros((ny, nx))
        self.p_star_new = np.zeros_like(self.p_star_old)
        self.p_prime = np.zeros_like(self.p_star_old)

        self.dx = mesh.dx
        self.dy = mesh.dy
        self.alpha_u = 0.7
        self.alpha_v = 0.3
        self.alpha_p = 0.3

    def calculate_diffusion_terms(self):
        De = Dw = self.mu / np.mean(self.dx)
        Dn = Ds = self.mu / np.mean(self.dy) if len(self.mesh.shape) == 2 else 0.0
        return De, Dw, Dn, Ds
    
    def generate_initial_guesses(self):
        is_2d = len(self.mesh.shape) == 2

        p_left = self.boundary.left_boundary['value'] if self.boundary.left_boundary and self.boundary.left_boundary['var'] == 'p' else None
        p_right = self.boundary.right_boundary['value'] if self.boundary.right_boundary and self.boundary.right_boundary['var'] == 'p' else None

        if p_left is not None and p_right is not None:
            xmin = self.mesh.domain.x_min
            xmax = self.mesh.domain.x_max
            for idx in self.mesh.interior_node_indices:
                if is_2d:
                    j, i = idx
                    node = self.mesh.nodes[j, i]
                else:
                    i = idx
                    j = 0
                    node = self.mesh.nodes[i]

                if node is None:
                    continue

                x = node.position[0]
                p_init = np.interp(x, [xmin, xmax], [p_left, p_right])
                self.p_star_old[j, i] = self.p_star_new[j, i] = p_init

        u_val = None
        if self.boundary.left_boundary and self.boundary.left_boundary['var'] == 'u':
            u_val = self.boundary.left_boundary['value']
        elif self.boundary.right_boundary and self.boundary.right_boundary['var'] == 'u':
            u_val = self.boundary.right_boundary['value']

        if u_val is not None:
            self.u_star_old.fill(u_val)
            self.u_star_new.fill(u_val)

        if is_2d:
            v_val = None
            if self.boundary.lower_boundary and self.boundary.lower_boundary['var'] == 'v':
                v_val = self.boundary.lower_boundary['value']
            elif self.boundary.upper_boundary and self.boundary.upper_boundary['var'] == 'v':
                v_val = self.boundary.upper_boundary['value']

            if v_val is not None:
                self.v_star_old.fill(v_val)
                self.v_star_new.fill(v_val)

        self.p_prime.fill(0.0)

                
    
    def solve_discretized_x_momentum(self):
        dx, dy, rho, mu = self.dx, self.dy, self.rho, self.mu
        De, Dw, Dn, Ds = self.calculate_diffusion_terms()
        is_2d = len(self.mesh.shape) == 2

        for j, i in self.mesh.u_face_indices:
            face = self.mesh.u_faces[j, i]

            # convective fluxes
            Fe = 0.5 * rho * (self.u_star_old[j, i+1] + self.u_star_old[j, i]) if i+1 < self.u_star_old.shape[1] else 0.0
            Fw = 0.5 * rho * (self.u_star_old[j, i] + self.u_star_old[j, i-1]) if i-1 >= 0 else 0.0
            Fn = 0.5 * rho * (self.v_star_old[j+1, i] + self.v_star_old[j+1, i-1]) if is_2d and j+1 < self.v_star_old.shape[0] else 0.0
            if is_2d and j < self.v_star_old.shape[0] - 1:  # Ensure j+1 is within bounds
             Fs = 0.5 * rho * (self.v_star_old[j+1, i] + self.v_star_old[j+1, i]) if i > 0 else 0.0
            else:
                Fs = 0.0
            # pressure gradient (pW - pE)
            pe = self.p_star_old[j, i] 
            pw = self.p_star_old[j, i-1] 

            face.u_old = self.u_star_old[j, i]
            face.calculate_x_coefficients(dx[i], dy[j], De, Fe, Dw, Fw, Dn, Fn, Ds, Fs, pe, pw, self.alpha_u, is_2d)

        # now solve A * u = b
        self.mesh.alpha_u = self.alpha_u
        self.mesh.alpha_v = self.alpha_v
        A, b = self.mesh.build_A_matrix('u')
        from scipy.sparse.linalg import spsolve
        u_flat = spsolve(A, b)
        for k, (j, i) in enumerate(self.mesh.u_face_indices):
            self.u_star_new[j, i] = u_flat[k]

            

        
    
    def solve_discretized_y_momentum(self):
        dx, dy, rho, mu = self.dx, self.dy, self.rho, self.mu
        De, Dw, Dn, Ds = self.calculate_diffusion_terms()
        is_2d = len(self.mesh.shape) == 2

        if not is_2d:
            return  # no v-momentum solve in 1D

        for j, i in self.mesh.v_face_indices:
            face = self.mesh.v_faces[j, i]

            # convective fluxes
            Fe = 0.5 * rho * (self.u_star_new[j, i+1] + self.u_star_new[j-1, i+1]) if i+1 < self.u_star_new.shape[1] else 0.0
            Fw = 0.5 * rho * (self.u_star_new[j, i] + self.u_star_new[j-1, i]) if i > 0 else 0.0
            Fn = 0.5 * rho * (self.v_star_old[j+1, i] + self.v_star_old[j, i]) if j+1 < self.v_star_old.shape[0] else 0.0
            Fs = 0.5 * rho * (self.v_star_old[j, i] + self.v_star_old[j-1, i]) if j > 0 else 0.0

            # pressure gradient (pN - pS)
            pn = self.p_star_old[j+1, i] if j+1 < self.p_star_old.shape[0] else 0.0
            ps = self.p_star_old[j, i] if j < self.p_star_old.shape[0] else 0.0

            face.v_old = self.v_star_old[j, i]
            face.calculate_y_coefficients(dx[i], dy[j], De, Fe, Dw, Fw, Dn, Fn, Ds, Fs, ps, pn, self.alpha_v, is_2d)

        # now solve A * v = b
        self.mesh.alpha_u = self.alpha_u
        self.mesh.alpha_v = self.alpha_v
        A, b = self.mesh.build_A_matrix('v')
        from scipy.sparse.linalg import spsolve
        v_flat = spsolve(A, b)
        for k, (j, i) in enumerate(self.mesh.v_face_indices):
            self.v_star_new[j, i] = v_flat[k]



    def solve_pressure_correction(self):
        dx, dy, rho = self.dx, self.dy, self.rho
        is_2d = len(self.mesh.shape) == 2

        for idx in self.mesh.interior_node_indices:
            if is_2d:
                j, i = idx
                node = self.mesh.nodes[j, i]
            else:
                i = idx
                j = 0
                node = self.mesh.nodes[i]

            if node is None:
                continue

            # Get face coefficients
            uE = self.mesh.u_faces[j, i+1] if i+1 < self.mesh.u_faces.shape[1] else None
            uW = self.mesh.u_faces[j, i] if i < self.mesh.u_faces.shape[1] else None
            vN = self.mesh.v_faces[j+1, i] if is_2d and j+1 < self.mesh.v_faces.shape[0] else None
            vS = self.mesh.v_faces[j, i] if is_2d and j < self.mesh.v_faces.shape[0] else None

            aP_u_E = uE.aP if uE is not None else 1e-12
            aP_u_W = uW.aP if uW is not None else 1e-12
            aP_v_N = vN.aP if vN is not None else 1e-12
            aP_v_S = vS.aP if vS is not None else 1e-12

            u_star_W = self.u_star_new[j, i] if i < self.u_star_new.shape[1] else 0.0
            u_star_E = self.u_star_new[j, i+1] if i+1 < self.u_star_new.shape[1] else 0.0
            v_star_S = self.v_star_new[j, i] if is_2d and j < self.v_star_new.shape[0] else 0.0
            v_star_N = self.v_star_new[j+1, i] if is_2d and j+1 < self.v_star_new.shape[0] else 0.0

            # Pass indexed values for dx and dy
            node.define_pressure_correction_coefficients(
                self.alpha_u, self.alpha_v, rho, dx[i], dy[j],
                aP_u_E, aP_u_W,
                aP_v_N, aP_v_S,
                u_star_W, u_star_E,
                v_star_S, v_star_N
            )

        A, b = self.mesh.build_A_matrix('p')
        from scipy.sparse.linalg import spsolve
        p_prime_flat = spsolve(A, b)

        for k, idx in enumerate(self.mesh.interior_node_indices):
            if is_2d:
                j, i = idx
            else:
                i = idx
                j = 0
            self.p_prime[j, i] = p_prime_flat[k]





    def calculate_new_pressure_field(self):
        for idx in self.mesh.interior_node_indices:
            # Handle 1D or 2D indexing properly
            if isinstance(idx, tuple):  # 2D case
                j, i = idx
            else:  # 1D case
                i = idx
                j = 0  # or set to 0 for 1D, or another logic if needed
            self.p_star_new[j, i] = self.p_star_old[j, i] + self.alpha_p * self.p_prime[j, i]

    def calculate_new_velocity_field(self):
        dx, dy, rho = self.dx, self.dy, self.rho
        is_2d = len(self.mesh.shape) == 2

        # Correct u
        for j, i in self.mesh.u_face_indices:
            if i+1 < self.p_prime.shape[1]:
                dp = self.p_prime[j, i-1] - self.p_prime[j, i]
            else:
                dp = 0.0

            aP = self.mesh.u_faces[j, i].aP if self.mesh.u_faces[j, i].aP else 1e-12

            # Use dy[j] to index the appropriate value of dy for each row j
            self.u_star_new[j, i] += dy[j] * dp / (rho * aP)

        # Correct v (only for 2D)
        if is_2d:
            for j, i in self.mesh.v_face_indices:
                if j+1 < self.p_prime.shape[0]:
                    dp = self.p_prime[j-1, i] - self.p_prime[j, i]
                else:
                    dp = 0.0

                aP = self.mesh.v_faces[j, i].aP if self.mesh.v_faces[j, i].aP else 1e-12

                # Use dx[i] to index the appropriate value of dx for each column i
                self.v_star_new[j, i] += dx[i] * dp / (rho * aP)



    def set_old_equal_to_new(self):
        self.u_star_old = self.u_star_new.copy()
        if self.v_star_old is not None:
            self.v_star_old = self.v_star_new.copy()
        self.p_star_old = self.p_star_new.copy()
    
    def print_and_plot_fields(self):
        # Check for 1D or 2D mesh
        is_2d = len(self.mesh.shape) == 2

        # --- u-face velocities ---
        print("--- u-face velocities ---")
        for j, i in self.mesh.u_face_indices:
            face = self.mesh.u_faces[j, i]
            if face is not None:
                print(f"u[{j}, {i}] = {face.u_old}")

        # --- pressure nodes ---
        print("--- pressure nodes ---")
        if is_2d:
            for j, i in self.mesh.interior_node_indices:
                p = self.p_star_new[j, i]  # For 2D mesh, use both indices
                x_position = self.mesh.nodes[j, i].position[0]  # Get the x-position of the node
                print(f"p[{j}, {i}] at x={x_position}: {p}")
        else:
            for i in range(self.mesh.shape[0]):  # For 1D mesh, iterate over the single dimension
                p = self.p_star_new[i]  # Access pressure for 1D
                x_position = self.mesh.nodes[i].position[0]  # Get the x-position of the node
                print(f"p[{i}] at x={x_position}: {p}")

        # Optionally, plot the fields here if needed








    def run(self):
        is_2d = len(self.mesh.shape) == 2
        
        self.generate_initial_guesses()
        converged = False
        max_iters = 2
        self.iterations = 0

        while not converged and self.iterations < max_iters:
            self.solve_discretized_x_momentum()
            if is_2d:
                self.solve_discretized_y_momentum()
            self.solve_pressure_correction()
            self.calculate_new_pressure_field()
            self.calculate_new_velocity_field()
            self.set_old_equal_to_new()

            # Convergence check
            residual = np.linalg.norm(self.p_prime)
            if residual < self.tolerance:
                converged = True

            self.iterations += 1
        self.print_and_plot_fields()

    