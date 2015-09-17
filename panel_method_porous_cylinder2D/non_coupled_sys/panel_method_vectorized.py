import numpy as np

def panel_method_vectorized(a, b, N, Z, tau=0.5, A = 1.0, rho=1.0):
    
    cz = 0.5 * (Z[1:] + Z[:-1])           # complex midpoint
    L = abs(Z[1:] - Z[:-1])               # Panel length
    n1 = (-Z.imag[1:] + Z.imag[:-1]) / L  # x-comp normal vector
    n2 = ( Z.real[1:] - Z.real[:-1]) / L  # y-comp normal vector
    
    # Matrix coefficients P and q
    P = np.zeros((N, N), dtype=complex)
    q = np.zeros((N, N), dtype=complex)
    
    # Opening angle given by the law of cosines
    for i in xrange(N):
        b = abs(cz[i] - Z[1:])
        c = abs(cz[i] - Z[:-1])
        q[i] = 0.5 * L * ( np.log(b) + np.log(c) ) # RHS integral
        P[i] = -np.arccos((b**2 + c**2 - L**2) / (2*b*c))
    np.fill_diagonal(P, -np.pi) # phi = -pi, for i = j
    
    # Parameters for fluid and porous media
    mu = 0.4                      # Discharge coefficient (Molin p. 3)
    K = (1 - tau) / (mu*tau**2)   # Resistance coefficient (Molin p. 3)
    KC = A / (2*a)                # Keulegan-Carpenter number (A / (2*a) )
    w = 1.4                       # Frequency omega
    
    # Additional term due to porosity
    P += 1j * w * q / (K * KC)
    
    # rhs
    Q = np.transpose([np.dot(q, n1), np.dot(q, n2), \
                      np.dot(q, (cz.real*n2 - cz.imag*n1))])

    # Velocity potential for each panel
    phi_i = np.linalg.solve(P, Q)

    # Added mass
    #d_phi_dn = 1j * w * d * n1 + rho - 1j * w * b / mu * Psi
    
    m11 = np.sum( phi_i.real[:, 0] * n1 * L )
    m22 = np.sum( phi_i.real[:, 1] * n2 * L )
    m12 = np.sum( phi_i.real[:, 0] * n2 * L )
    m66 = np.sum( phi_i.real[:, 2] * (cz.real*n2 - cz.imag*n1) * L )

    # Damping coefficients
    b11 = np.sum( -phi_i.imag[:, 0] * n1 * L / w)

    # Added mass and damping
    d_phi_dn = 1j * w * n1# - rho * 1j * w * b / mu * phi_i[:,0]
    f11 = -rho * np.sum(phi_i[:,0] * d_phi_dn * L)

    m11 = -f11.real/w**2
    b11 = -f11.imag/w
    
    return m11, m22, m12, m66, b11
