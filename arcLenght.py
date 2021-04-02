"""
Program that returns the arc length from 'a' to 'b', of the parametrized function alpha(t)

alpha(t) function is parametrized, but not neccesary arc's length parametrized
"""

# Libraries
import sympy as sp
import sympy.plotting as spp

# Symbols to be used, do not change these
t = sp.symbols('t', positive=True)

# Function alpha(t)=([[alpha1(t), alpha2(t), alpha3(t)]])
alpha = sp.Matrix([[sp.cos(t), sp.sin(t), t]])

# Tangent vector (d/dt)alpha(t)
da = sp.simplify(sp.diff(alpha, t))
print("alpha(t)'s Tangent vector T(t)=(d/dt)alpha(t): T(t)=%s \n" %
      (da))

# Norm ||(d/dt)alpha(t)||
nda = sp.simplify(sp.sqrt(da.dot(da)))
print("Norm of alpha(t): ||(d/dt)alpha(t)||=%s \n" % nda)

# Arc Lenght. Int_{t=a}^{b} alpha(t) dt
# Modify this to change the limits of integration
a = -3  # Inferior limit
b = 1  # Superior limit
try:
    """
    This try-except checks if the integral can be solved
    - If can be solved returns the arc length from 'a' to 'b'
    - Otherwise, prints 'can not be solved
    """
    Lab = sp.simplify(sp.integrate(nda, (t, a, b)))
    print("ArcLenght of alpha(t): %s" % (Lab))
except:
    print("The integral can not be solved")


# Plotting section
# The domain of the plots is the same of the limits of integration in 'Lab', you can change these if you want

# Obtaining the components of alpha(t)
x = alpha[0]
y = alpha[1]
z = alpha[2]

# 2d
# spp.plot_parametric(x,y,(t,a,b),title=r"Plot of $\alpha(t)$")

# 3d
spp.plot3d_parametric_line(x, y, z, (t, a, b), title=r"Plot of $\alpha(t)$")
