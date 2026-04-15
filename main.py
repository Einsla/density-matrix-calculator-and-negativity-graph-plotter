import negativity_plotter as neg
import density_calculation as dc
import sympy as sp

pi = sp.pi
sqrt = sp.sqrt

# initial_state = 0.5 * sp.Matrix([[1.0, 0.0, 0.0, 1.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 1.0]])
a = [0.5, 0.5, 0.5, 0.5]
v = [pi/2, pi, -pi, -pi/2]
initial_state = dc.ket_to_density(a, v)
# initial_state = sp.Matrix([
#     [0.9975, 0.0, 0.0, 0.0499],
#     [0.0, 0.0, 0.0, 0.0],
#     [0.0, 0.0, 0.0, 0.0],
#     [0.0499, 0.0, 0.0, 0.0025]])
neg.plot_negativity("continuous", 80, initial_state, 0.1, 'MachinePrecision', 10)
dc.close()