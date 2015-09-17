import numpy as np, sys

def panel_method_vectorized(a, b, N, Z, tau=0.5, d=1.0, rho=1.0):

    # d   - amplitude of oscillation
    # tau - porosity parameter, 0 is solid, 1 is open
    
    cz = 0.5 * (Z[1:] + Z[:-1])           # complex midpoint
    L = abs(Z[1:] - Z[:-1])               # Panel length
    n1 = (-Z.imag[1:] + Z.imag[:-1]) / L  # x-comp normal vector
    n2 = ( Z.real[1:] - Z.real[:-1]) / L  # y-comp normal vector
    
    #
    theta = np.zeros( (N, N) )
    intLog = np.zeros( (N, N) )
    A = np.pi * np.eye(2*N, dtype=complex)
    
    for i in xrange(N):
        b = abs(cz[i] - Z[1:])
        c = abs(cz[i] - Z[:-1])
        intLog[i] = 0.5 * L * ( np.log(b) + np.log(c) )
        # Opening angle by the law of cosines:
        theta[i] = np.arccos((b**2 + c**2 - L**2) / (2*b*c))
    np.fill_diagonal(theta, np.pi)

    # Parameters for fluid and porous media
    mu = 0.4                      # Discharge coefficient (Molin p. 3)
    K = (1 - tau) / tau**2        # Resistance coefficient (Molin p. 3)
    KC = d / (2*a)                # Keulegan-Carpenter number
    b = 1 / (K * KC)              # Porosity coefficient from Chwang
    w = 1.0                       # Frequency omega
    
    #
    A[N:, :N] = theta
    A[:N, N:] = theta
    A[N:, N:] += 2 * 1j * rho * w * b / mu * intLog

    # Rhs
    B = np.zeros(2*N, dtype=complex)
    B[N:] = -2 * 1j * w * np.dot(intLog, n1) # * d

    # Psi and Phi on each panel
    PHI = np.linalg.solve(A, B)
    Phi = PHI[:N]
    Psi = PHI[N:]
    
    # Inner and outer potential
    phi_e = 0.5 * (Phi + Psi)
    phi_i = 0.5 * (Phi - Psi)

    # Added mass and damping
    d_phi_dn = 1j * w * n1 - rho * 1j * w * b / mu * Psi
    f11 = -rho * np.sum(Phi * d_phi_dn * L)

    # Drag coefficient
    U = 1                               # Stream with unit speed
    S = np.pi * a                       # Area of cylinder
    Fd = -rho * 1j * w * np.sum(Phi*n1) # Total drag force
    Cd = 2 * Fd / (rho * U**2 * S)      # Drag coeffient
    
    return f11.real / w**2, -f11.imag / w, Cd
