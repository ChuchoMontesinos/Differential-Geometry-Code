"""
Program that receives the alpha(t) curve, checks if it is parametrized by arc lenght
- If it is parametrized by arc length approximates the function alpha(s) with the Frenet-Serret Equations using Taylor series following the formula:
    alpha(s)=alpha(0)+s*Tangent(s)+((s**2)/2))*Curvature(s)*Normal(s)+((s**3)/6)*Curvature(s)*Torsion(s)*Binormal(s)

-Otherwise, reparametrizes the curve and the approximates the function
"""

# Importing the libraries
import sympy as sp

# Declaring the symbols
t, s, u = sp.symbols('t s u', positive=True)

# alpha(t) function,
alpha = sp.Matrix([[sp.sin(t), sp.cos(t), t]])
print("alpha(t)=%s" % (alpha))

# derivative of alpha(t)
da = sp.simplify(alpha.diff(t))
print("(d/dt)alpha(t)=%s" % (da))

# Norm of the derivative
nda = sp.simplify(sp.sqrt(da.dot(da)))
print("||(d/dt)alpha(t)||=%s" % (nda))


# Checks if alpha(t) is parametrized by arc's lenght
if nda == 1:
    print("alpha(t) is parametrized by arc's lenght")
    aps = alpha.subs(t, s)  # substituting t->s
    S = aps.subs(s, 0)  # Evaluating aplha(s) when s=0
    for b in range(1, 4):
        # This 'for' computes the 3 first terms of the Taylor series
        # expression of the taylor expansion
        da = (1/sp.factorial(b))*(s**b)*sp.diff(aps, (s, b))
        # Printing the expressions
        print("Derivative %s: %s" % (b, sp.simplify(da)))
        S += da  # Adding the expressions
    print("Adding the results we get: %s\n" % (sp.simplify(S)))
    print("Substituting t=0: %s" % (S.subs(s, 0)))
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
            print("t=:%s\n" % (desp))
            # Reparametrization of alpha(t)
            nalpha = sp.simplify(alpha.subs(t, desp[0]))
            print("alpha(s): %s" % (nalpha))
            Sp = nalpha.subs(s, 0)  # Initializing the initial value
            for b in range(1, 4):
                # This 'for' computes the 3 first terms of the Taylor series
                # expression of the taylor expansion
                da = (1/sp.factorial(b))*(s**b)*sp.diff(nalpha, (s, b))
                # Printing the expressions
                print("Derivative %s: %s" % (b, sp.simplify(da)))
                Sp += da  # Adding the expressions
            print("Adding the results we get: %s\n" % (sp.simplify(Sp)))
            print("Substituting t=0: %s" % (Sp.subs(s, 0)))

        except:
            print("The arc's length function can not be solved")
            print("alpha(s) can not be parametrized by arc's length")

    except:
        # This code is executed if 'lenArc' could not be solved
        print("The integral can not be solved\n")
        print("alpha(t) can not be parametrized by arc's length")
