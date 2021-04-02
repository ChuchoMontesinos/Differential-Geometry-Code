"""
Program that plots the function alpha(t)
According to alpha(t) you can uncomment the code of the plotting section
"""

# Libraries
import sympy as sp
import sympy.plotting as spp

# Symbols to be used, do not change these
t = sp.symbols('t', positive=True)

# Constants

# alpha(t) vector
alpha = sp.Matrix([[t, t**2]])


# Plotting section

# Obtainting the components of alpha(t)
x = alpha[0]
y = alpha[1]
# z=alpha[2]

# 2D
spp.plot_parametric(x, y, (t, 0, 1), title=r"Plot of $\alpha(t)$")

# 3D
#spp.plot3d_parametric_line(x,y,z,(t,-3,3),title=r"Plot of $\alpha(t)$")
