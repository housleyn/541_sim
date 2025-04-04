import matplotlib.pyplot as plt
import numpy as np

class Domain():
    def __init__(self):

        self.lower_boundary_func = None
        self.upper_boundary_func = None
        self.left_boundary_func = None
        self.right_boundary_func = None

    def define_lower_boundary(self, func):
        self.lower_boundary_func = func

    def define_upper_boundary(self, func):
        self.upper_boundary_func = func

    def define_left_boundary(self, func):
        self.left_boundary_func = func

    def define_right_boundary(self, func):
        self.right_boundary_func = func

    def is_inside(self, x, y):
        if self.lower_boundary_func(x, y) < 0: return False  # below lower boundary
        if self.upper_boundary_func(x, y) > 0: return False  # above upper boundary
        if self.left_boundary_func(y, x) > x: return False   # left of left
        if self.right_boundary_func(y, x) < x: return False  # right of right
        return True


    def plot_domain(self, x_range, y_range):
        x = np.linspace(x_range[0], x_range[1], 200)
        y_lower = np.array([self.lower_boundary_func(xi, 0) for xi in x])
        y_upper = np.array([self.upper_boundary_func(xi, 0) for xi in x])

        y_left_min = self.lower_boundary_func(x_range[0], 0)
        y_left_max = self.upper_boundary_func(x_range[0], 0)
        y_right_min = self.lower_boundary_func(x_range[1], 0)
        y_right_max = self.upper_boundary_func(x_range[1], 0)

        y_left = np.linspace(y_left_min, y_left_max, 200)
        y_right = np.linspace(y_right_min, y_right_max, 200)
        x_left = np.array([self.left_boundary_func(0, yi) for yi in y_left])
        x_right = np.array([self.right_boundary_func(0, yi) for yi in y_right])

        plt.plot(x, y_lower, 'k', label='Lower Boundary')
        plt.plot(x, y_upper, 'k', label='Upper Boundary')
        plt.plot(x_left, y_left, 'k', label='Left Boundary')
        plt.plot(x_right, y_right, 'k', label='Right Boundary')

        plt.xlim(x_range)
        plt.ylim(y_range)
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')
        plt.title('Domain with Boundaries')
        plt.grid()
        plt.axis('equal')
        plt.show()




