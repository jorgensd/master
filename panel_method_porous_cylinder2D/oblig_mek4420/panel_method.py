import numpy as np, loop # Migrated fortran loop

def panel_method(a, b, N, Z):

    cz = 0.5 * (Z[1:] + Z[:-1])           # complex midpoint
    L = abs(Z[1:] - Z[:-1])               # Panel length
    n1 = (-Z.imag[1:] + Z.imag[:-1]) / L  # x-comp normal vector
    n2 = ( Z.real[1:] - Z.real[:-1]) / L  # y-comp normal vector
        
    # P*phi = Q
    P = np.zeros( (N, N), order='Fortran' )
    Q = np.zeros( (N, 3), order='Fortran' )
    
    # Migrated fortran loop computes the lhs matrix P,
    # and the rhs integral for surge, heave and pitch.
    P, Q = loop.loop(N, Z, cz, L, n1, n2, P, Q)
        
    # Solve the linear matrix system, P * phi = Q, yields phi on each panel
    phi_i = np.linalg.solve(P, Q)
    
    # Added mass.
    m11 = np.sum( phi_i[:, 0] * n1 * L )
    m22 = np.sum( phi_i[:, 1] * n2 * L )
    m12 = np.sum( phi_i[:, 0] * n2 * L )
    m66 = np.sum( phi_i[:, 2] * (cz.real*n2 - cz.imag*n1) * L )
        
    return m11, m22, m12, m66
