import numpy as np, sys
import cmath as cm

def panel_method_vectorized(a, b, N, Z):
    
    cz = 0.5 * (Z[1:] + Z[:-1])           # complex midpoint
    cx = cz.real                          # x-coor panel midpoint
    cy = cz.imag                          # y-coor panel midpoint
    L = abs(Z[1:] - Z[:-1])               # Panel length
    n1 = (-Z.imag[1:] + Z.imag[:-1]) / L  # x-comp normal vector
    n2 = ( Z.real[1:] - Z.real[:-1]) / L  # y-comp normal vector
    
    # Matrix coefficients P and q
    P = np.zeros([N, N])
    q = np.zeros([N, N])

    # Opening angle given by the law of cosines
    for i in xrange(N):
        b = abs(cz[i] - Z[1:])
        c = abs(cz[i] - Z[:-1])
        P[i] = -np.arccos((b**2 + c**2 - L**2) / (2*b*c))
        q[i] = 0.5 * L * ( np.log(b) + np.log(c) ) # RHS integral
        #P[i] = np.angle(Z[1:] - cz[i]) - np.angle(Z[:-1] - cz[i])  

    P[np.isnan(P)] = 0 # Needed for rectangle
    np.fill_diagonal(P, -np.pi) # phi = -pi, for i = j
    
    # RHS
    Q = np.transpose([np.dot(q,n1),np.dot(q,n2),np.dot(q,(cx*n2 - cy*n1))])

    # Velocity potential for each panel
    phi_i = np.linalg.solve(P, Q)
    print phi_i[:, 0]
    
    # Added mass
    m11 = np.sum( phi_i[:, 0] * n1 * L )       # m11 for the whole body
    m22 = np.sum( phi_i[:, 1] * n2 * L )
    m12 = np.sum( phi_i[:, 0] * n2 * L )       # Cross coupling
    m66 = np.sum( phi_i[:, 2] * (cx*n2 - cy*n1) * L )

    return m11, m22, m12, m66
