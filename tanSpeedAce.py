"""
Program that returns the tangent vector, speed and acceleration of the parametrized function alpha(t) and plots alpha(t)

alpha(t) is parametrized, but not necessarily parameterized by arc length
"""

# Importing the libraries
import sympy as sp
import sympy.plotting as spp

# Symbols to be used, do not change this
t = sp.symbols('t', positive=True)

# alpha(t) function alpha(t)=([[alpha1(t), alpha2(t), alpha3(t)]])
alpha = sp.Matrix([[sp.cos(t), sp.sin(t), sp.exp(t)]])

# Tangent vector (d/dt)alpha(t)
da = sp.simplify(sp.diff(alpha, t))
print("alpha(s)'s Tangent vector T(t)=(d/dt)alpha(t): T(t)=%s \n" %
      (da))

# Speed ||(d/dt)alpha(t)||
print("Speed of alpha(t): ||(d/dt)alpha(t)||=%s\n" %
      (sp.simplify(sp.sqrt(da.dot(da)))))

# Acceleration ((d2)/d(t2))alpha(t)
dda = sp.simplify(da.diff(t))
print("Acceleration of alpha(t) ((d2)/(dt2))alpha(t): a(t)=%s" %
      (dda))

# Plotting

"""
Plots alpha(t) uncomment the option according to alpha(t)
The code inside $'s signs is LaTeX code
"""

x = alpha[0]
y = alpha[1]
z = alpha[2]

# 2d
#spp.plot_parametric(x, y, (t, -2, 2), title=r"Plot of $\alpha(t)$")

# 3d
spp.plot3d_parametric_line(x, y, z, (t, -2, 2), title=r"Plot of $\alpha(t)$")
