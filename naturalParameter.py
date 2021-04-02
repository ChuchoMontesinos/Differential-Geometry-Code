"""
Program that returns the natural parametrization of the function alpha(t). 
This takes the alpha(t) function and returns the alpha(s) where 's' is obteined by the arc's length function 
"""

# Libraries
import sympy as sp
import sympy.plotting as spp

"""
Symbols to be used during the process, you must not change this
't' is the original parameter
'u' is the variable to reparametrize alpha(t) to solve the arc's length function
's' denotes the arc's lenght function
"""
t, u, s = sp.symbols('t u s', positive=True)

# Auxiliar symbols

# Constants that appear in the alpha(t) function, you can change this according to your problem
# If you change these nothing is going to be affected


# alpha(t) function
alpha = sp.Matrix([[t, t**2]])

# alpha(t)'s derivative (d/dt)alpha(t)
da = sp.simplify(sp.diff(alpha, t))
print("Derivative of alpha(t): (d/dt)alpha(t)=%s \n" %
      (da))

# Norm of the derivative ||(d/dt)alpha(t)||
nda = sp.simplify(sp.sqrt(da.dot(da)))
print("Norm of Derivative: %s\n\n" % (nda))

# Process to verify if alpha(t) is parametrized by arc's length

if nda == 1:
    print("alpha(t) is parametrized by arc length, so, alpha(t)=alpha(s)")
else:
    print("alpha(t) is not parametrized by arc lenght\n")

    try:
        """
        This try-except is to verify if the integral can be solved, 
        - If the integral can be solved then the next step is get the expression explicit in 's'
        - Otherwise, prints 'can not be solved'
        """
        # arc's length function, Int_{u=0}^{t} alpha(u) du
        lenArc = sp.simplify(sp.integrate(nda.subs(t, u), (u, 0, t)))
        print("Arc's lenght function s(t): %s\n" % (lenArc))

        # If the expression above could be solved, the next step is try to get the implicit function in 's'
        equa = sp.Eq(s, lenArc)
        print("Equation to be solved: %s\n " % (equa))
        try:
            """
            This try-except checks if the expression in 'equa' can be solved
            - If the expression can be solved reparametrices alpha(t) by arc's lenght
            - Otherwise, prints 'can not be parametrized'
            """
            # Print the function explicit in 's'
            desp = sp.solve(equa, t)
            print("Parameter f(s):%s\n" % (desp))
            # Reparametrization of alpha(t)
            nalpha = alpha.subs(t, desp[0])
            print("alpha(s): %s" % (nalpha))

        except:
            print("The arc's length function can not be solved")
            print("alpha(s) can not be parametrized by arc's length")

    except:
        # This code is executed if 'lenArc' could not be solved
        print("The integral can not be solved\n")
        print("alpha(t) can not be parametrized by arc's length")


# Plotting section

x = alpha[0]
y = alpha[1]
#z = alpha[2]

# 2d
spp.plot_parametric(x, y, (t, 0, 1), legend=True)

# 3d
#spp.plot3d_parametric_line(x,y,z,(t,-2*sp.pi,2*sp.pi), label=True)
