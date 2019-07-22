from sympy import *

# Primitive variables
rho, u, e, gamma = symbols('rho u e gamma',real=True)

# Conserved variables
q1, q2, q3 = symbols('q1 q2 q3')

# Vector of primitive variables
P = Matrix([rho, u, e])

# Vector of conserved variables
E = e + (1/2)*u**2 # Total energy
Q = Matrix([rho, rho*u, rho*E])

# Flux vector
p = (gamma - 1) * (rho*e) # Ideal gas law
F = Matrix([rho*u, rho*u**2 + p, rho*u*E + p*u])

# Compute flux Jacobian from chain rule:
# dF/dQ = (dF/dP)(dP/dQ) = (dF/dP)(dQ/dP)^-1
dFdP = F.jacobian(P)
dQdP = Q.jacobian(P)
dFdQ = dFdP * dQdP**(-1)

# Rewrite using conserved variables
dFdQ = dFdQ.subs(e, q3/rho - (1/2)*u**2)
dFdQ = dFdQ.subs(rho, q1)
dFdQ = dFdQ.subs(u, q2/q1)
dFdQ.simplify()

# Diagonalize
#T, Lambda = dFdQ.diagonalize()
#T_inv = T**(-1)
T = ''
T_inv = ''
Lambda = simplify(dFdQ.eigenvals())

# Output results
print("""

Primitive Variables:

""")
print(latex(P))

print("""

Conservative Variables:

""")
print(latex(Q))

print("""

Flux Vector:

""")
print(latex(F))

print("""

Jacobian of Flux Vector w.r.t. Conserved Variables:

""")
print(latex(dFdQ))

print("""

Eigenvalues of the Jacobian of the Flux Vector:

""")
print(latex(Lambda))

print("""

Right Eigenvalues of the Jacobian of the Flux Vector:

""")
print(latex(T))

print("""

Left Eigenvalues of the Jacobian of the Flux Vector:

""")
print(latex(T_inv))

print()
