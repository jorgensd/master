import numpy as np, loop # Migrated fortran loop

def panel_method(a, b, N, Z, tau=0.5, d=1, rho=1, w=1):

    #
    cz = 0.5 * (Z[1:] + Z[:-1])           # complex midpoint
    L = abs(Z[1:] - Z[:-1])               # Panel length
    n1 = (-Z.imag[1:] + Z.imag[:-1]) / L  # x-comp normal vector
    n2 = ( Z.real[1:] - Z.real[:-1]) / L  # y-comp normal vector

    # Parameters for fluid and porous media
    mu = 0.4                      # Discharge coefficient (Molin p. 3)
    K = (1 - tau) / (mu*tau**2)   # Resistance coefficient (Molin p. 3)
    KC = d / (2.*a)                # Keulegan-Carpenter number
    b = 1 / (K * KC)              # Porosity coefficient from Chwang
    w = 1.0                       # Frequency omega
    
    #
    theta = np.zeros((N, N), dtype=complex, order='Fortran')
    intLog = np.zeros((N, N), dtype=complex, order='Fortran')
    A = np.pi * np.eye(2*N, dtype=complex)
    B = np.zeros(2*N, dtype=complex, order='Fortran')
                
    #
    B = loop.loop(N, Z, cz, L, n1, theta, intLog, A, B, mu, tau, a, d, w, rho)
            
    # Psi and Phi on each panel
    Phi = B[:N]
    Psi = B[N:]
    
    # Inner and outer potential
    phi_e = 0.5 * (Phi + Psi)
    phi_i = 0.5 * (Phi - Psi)

    # Added mass
    m11 = np.sum(phi_e.real * n1 * L) + np.sum(phi_i.real * n1 * L)
    
    # Damping coefficient
    b11 = np.sum(phi_e.imag * n1 * L / w) + np.sum(phi_i.imag * n1 * L / w)

    return m11, b11
