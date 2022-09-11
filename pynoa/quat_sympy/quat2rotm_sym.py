from sympy import zeros
def quat2rotm_sym(q):
    # Quaternion to rotation matrix converter. Also accepts symbolic input
    # qw = q[0] 
    # qx = q[1]
    # qy = q[2]
    # qz = q[3]
    qw = q[0] 
    qx = q[1]
    qy = q[2]
    qz = q[3]
    R = zeros(3)
    R[0,0] = qw**2 + qx**2 - qy**2 - qz**2
    R[1,0] = 2*qx*qy + 2*qw*qz
    R[2,0] = 2*qx*qz - 2*qw*qy
    R[0,1] = 2*qx*qy - 2*qw*qz
    R[1,1] = qw**2 - qx**2 + qy**2 - qz**2
    R[2,1] = 2*qy*qz + 2*qw*qx
    R[0,2] = 2*qx*qz + 2*qw*qy
    R[1,2] = 2*qy*qz - 2*qw*qx
    R[2,2] = qw**2 - qx**2 - qy**2 + qz**2
    return R
