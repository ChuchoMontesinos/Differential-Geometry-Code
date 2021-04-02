"""
Program that returns the Frenet-Serret equations. You introduce the alpha(t) curve parametrized
First checks if the curve is parametrized by arc's length
- If it is parametrized by arc's length calculates the Frenet-Serret equations,
- Otherwise reparametrizes the curve and then calculates the Frenet-Serret equations
"""

# Libraries
import sympy as sp

# Functions


def tangent(Q):
    """
    Function that returns the tangent vector T(s)=(d/ds)Q(s)
    """
    return sp.simplify(Q.diff(s))


def curvature(Q):
    """
    Function that returns the curvature kappa(s)=||(d/ds)Q(s)||
    """
    # Derivative of Q(s)
    T = Q.diff(s, 2)
    return sp.simplify(sp.sqrt(T.dot(T)))


def binormal(Tv, Nv):
    """
    Function that returns the binormal vector of Q(s), B(s)=Tv(s)xNv(s)
    'Tv(s)' is the tangent vector and 'Nv(s)' the normal vector
    """
    return sp.simplify(Tv.cross(Nv))


def torsion(Q):
    """
    Function that returns the torsion tau(s) of Q(s)
    """
    # Derivatives of alpha(s)
    dQ = Q.diff(s)
    ddQ = Q.diff(s, 2)
    dddQ = Q.diff(s, 3)
    # Cros product dQxddQ
    pc = dQ.cross(ddQ)
    # Triple escalar product
    tpe = sp.simplify(dQ.dot(ddQ.cross(dddQ)))
    npc = sp.simplify(pc.dot(pc))
    return sp.simplify(tpe/npc)


# Symbols to be used during the process. Do no change these
t, s, u = sp.symbols('t s u', positive=True)

# alpha(t) Function
alpha = sp.Matrix([[sp.cos(t), sp.sin(t), t]])

# Derivative (d/dt)alpha(t)
da = sp.simplify(sp.diff(alpha, t))
print("The derivative of alpha(t): (d/dt)alpha(t)=%s" % (da))

# Norm |(d/dt)alpha(t)|
nda = sp.simplify(sp.sqrt(da.dot(da)))
print("The norm of the derivative of alpha: ||(d/dt)alpha(t)||=%s" % (nda))

if nda == 1:
    print("alpha(t) is parametrized by arc's length, so alpha(t)=alpha(s)\n")

    # Changing the parameter alpha(t) -> alpha(s)
    alphas = alpha.subs(t, s)

    # Tangent vector (d/ds)alpha(s)=T(s)
    Ts = tangent(alphas)
    print("alpha(s)'s Tangent vector T(s)=(d/ds)alpha(s): T(s)=%s" %
          (Ts))

    # Derivative of Tangent vector (d/ds)T(s)
    dTs = sp.simplify(Ts.diff(s))
    print("Tangent's vector derivative (d/ds)T(s)=%s" %
          (dTs))

    # Curvature kappa(s)=|(d/ds)T(s)|
    kappas = curvature(alphas)
    print("alpha(s)'s curvature kappa(s)=|(d/ds)T(s)|: kappa(s)= %s " % (kappas))

    # Curvature ratio rho(s)=(1/kappa(s))
    print("alpha(s)'s ratio of curvature rho(s)=(1/kappa(s))): rho(s)= %s" %
          (1/kappas))

    # Normal vector N=(1/kappa(s))(d/ds)T -> N=(d/ds)T
    print("alpha(s)'s Normal vector N(s)=(1/kappa)((d/ds)T): N(s)=%s" %
          (dTs))

    # Binormal vector B=TxN
    Bs = binormal(Ts, dTs)
    print("alpha(s)'s Binormal vector B(s)=TxN: B(s)=%s" %
          (Bs))

    # Binormal vector derivative (d/ds)B(s)
    dBs = sp.simplify(Bs.diff(s))
    print("Binormal's vector derivative (d/ds)B(s): (d/ds)B(s)=%s" %
          (dBs))

    # Torsion tau(s)
    taus = torsion(alphas)
    print("alpha(s)'s torsion: tau(s)=%s" % (taus))

    if taus == 0:
        print("There is not torsion ratio sigma(s) because tau(s)=0")
    else:
        # Ratio of torsion: sigma(s)=1/tau
        print("Torsion's ratio sigma(s)=1/tau: sigma(s)=%s" % (1/taus))

else:
    print("alpha(t) is not para metrized by arc's length\n")

    # Changing the parameter alpha(s) -> alpha(u)
    au = alpha.subs(s, u)
    print("alpha(u): alpha(u)=%s" % (au))

    # derivative of alpha(u)
    dau = au.diff(u)
    print("The derivative of alpha(u): (d/du)alpha(s)=%s" % (dau))

    # Norm of the derivative |(d/du)alpha(u)|
    ndau = sp.simplify(sp.sqrt(dau.dot(dau)))
    print("alpha(u)'s derivative norm: |(d/du)alpha(u)|= %s" % (ndau))

    try:
        """
        This try-except tries to calculate the integral
        - If the integral can be calculated, the next step is to solve the equation
        - Otherwise, prints 'Integral can not be done'
        """
        # Arc's length function Int_{u=0}^{t} ||(d/du)alpha(u)|| du
        ALf = sp.simplify(sp.integrate(ndau, (u, 0, t)))
        print(
            "Arc's Length Function result s(t)=Int_{u=0}^{t} ||(d/du)alpha(u)|| du: s(t)=%s" % (ALf))

        # Obtaining the implicit expression for s(t)
        equa = sp.Eq(s, ALf)
        print("The equation to be solved is %s" % (equa))
        try:
            """
            This try-except tries to solve the Arc's Length function
            - If it can be solved then we can continue with the Frenet-Serret Equations
            - Otherwise prints "The function can not be reparametrized"
            """
            # Solving the equation
            sol = sp.solve(equa, t)
            print("Obtaining t(s): t=%s" % (sol[0]))

            # Replacing it on alpha(t)
            alphasr = sp.simplify(alpha.subs(t, sol[0]))
            print("The arc's lenght reparametrization is: alpha(s)=%s" %
                  (alphasr))

            # alpha(s)'s Tangent vector
            Tr = tangent(alphasr)
            print("alpha(s)'s Tangent vector (d/ds)alpha(s): T(s)=%s" %
                  (Tr))

            # Derivative of the tangent vector (d/ds)T(s)
            dTr = Tr.diff(s)
            print("alpha(s)'s Tangent vector derivative (d/ds)T(s): dT(s)=%s" %
                  (dTr))

            # alpha(s) curvature kappa
            kappar = curvature(alphasr)
            print("alpha(s)'s curvature kappa(s)= |(d/ds)T(s)|: kappa(s)=%s" % (kappar))

            # ratio of curvature rho(s)=1/kappa
            print("The curvature's ratio is: rho(s)=%s" % (1/kappar))

            # Normal vector N=(1/kappa(s))(d/ds)T -> N=(d/ds)T
            print("alpha(s)'s Normal vector N(s)=(1/kappa)((d/ds)T): N(s)=%s" %
                  (dTr))

            # Binormal vector B=TxN
            Br = binormal(Tr, dTr)
            print("alpha(s)'s Binormal vector B(s)=TxN: B(s)=%s" %
                  (Br))

            # Binormal vector derivative (d/ds)B(s)
            dBr = sp.simplify(Br.diff(s))
            print("Binormal's vector derivative (d/ds)B(s): (d/ds)B(s)=%s" %
                  (dBr))

            # Torsion tau(s)
            taur = torsion(alphasr)
            print("alpha(s)'s torsion: tau(s)=%s" % (taur))

            if taur == 0:
                print("There is not torsion ratio sigma(s) because tau(s)=0")
            else:
                # Ratio of torsion: sigma(s)=1/tau
                print("Torsion's ratio sigma(s)=1/tau: sigma(s)=%s" % (1/taur))

        except:
            print(
                "The equation can not be solved, alpha(t) can not be re parametrized by arc's length ")
    except:
        print(
            "Integral can not be solved, alpha(t) can not be parametrized by arc's length")
