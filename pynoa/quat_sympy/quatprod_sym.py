from sympy import Matrix
def quatprod_sym(q, r):
    # Custom Python script to calculate product of 2 quaternions
    q0=q[0] 
    q1=q[1] 
    q2=q[2] 
    q3=q[3]
    r0=r[0]
    r1=r[1] 
    r2=r[2] 
    r3=r[3]
    n0 = r0*q0 - r1*q1 - r2*q2 - r3*q3
    n1 = r0*q1 + r1*q0 - r2*q3 + r3*q2
    n2 = r0*q2 + r1*q3 + r2*q0 - r3*q1
    n3 = r0*q3 - r1*q2 + r2*q1 + r3*q0
    n = Matrix([[n0, n1, n2, n3]]).T
    return n