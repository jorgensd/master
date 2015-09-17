import numpy as np, sys

def panel_method_vectorized(a, b, N, Z, bpor=0.5, d=0.1,
                            rho=1.0, mu=0.3, w=0.5, ngc='64'):

    # d    - amplitude of oscillation
    # bpor - Porosity coefficient from Chwang
    # mu   - Discharge coefficient (Molin p. 3)
    # w    - Frequency omega

    cz = 0.5 * (Z[1:] + Z[:-1])           # complex midpoint
    L = abs(Z[1:] - Z[:-1])               # Panel length
    n1 = (-Z.imag[1:] + Z.imag[:-1]) / L  # x-comp normal vector
    n2 = ( Z.real[1:] - Z.real[:-1]) / L  # y-comp normal vector
    u0 = 1 #1j * w * d                       # complex velocity
    cx = cz.real
    cy = cz.imag
    x = Z.real
    y = Z.imag
    
    # Auxiliary matrices
    theta = np.zeros( (N, N) )
    intLog = np.zeros( (N, N) )
    A = np.pi * np.eye(2*N, dtype=complex)

    # Converts the rhs integral into a line integral, which can be
    # evaluated by the gaussian quadrature rule.
    def GQ(x0, x1, xm, y0, y1, ym, t):
        r = np.sqrt((0.5*(x0*(1-t)+x1*(1+t))-xm)**2+(0.5*(y0*(1-t)+y1*(1+t))-ym)**2)
        v = 0.5 * np.sqrt((x1 - x0)**2 + (y1 - y0)**2)
        return np.log(r)*v

    # Importing weight and points for gauss quad int 
    from readdata import readdata
    w_, p_ = readdata('weightsAndPoints%s.txt' % ngc)
    
    # Calculating the opening angle, theta, and solving the rhs integral
    # using 64 point Gauss method.
    for i in xrange(N):
        b = abs(cz[i] - Z[1:])
        c = abs(cz[i] - Z[:-1])
        theta[i] = np.arccos((b**2 + c**2 - L**2) / (2*b*c))
        #intLog[i] = 0.5 * L * ( np.log(b) + np.log(c) ) # Trapz
        for (wi, pi) in zip(w_, p_):
            intLog[i] += wi * GQ(x[:-1], x[1:], cx[i], y[:-1], y[1:], cy[i], pi)
    np.fill_diagonal(theta, np.pi)
    
    # See appendix A
    A[N:, :N] = theta    # upper right corner
    A[:N, N:] = -theta   # lower left corner
    A[N:, N:] -= 2 * 1j * rho * w * bpor / mu * intLog # lower right

    # Rhs
    B = np.zeros(2*N, dtype=complex)
    B[N:] = -2 * u0 * np.dot(intLog, n1)

    # Psi and Phi on each panel
    PHI = np.linalg.solve(A, B)
    Phi = PHI[:N]
    Psi = PHI[N:]

    # Inner and outer potential
    # phii = 0.5 * (Phi - Psi)
    # phie = 0.5 * (Phi + Psi)
    
    # Added mass and damping
    f11 = rho * np.sum(Psi * n1 * L) # dphidn = n1

    # Drag coefficient
    U = 1                               # Stream with unit speed
    S = np.pi * a                       # Area of cylinder
    Fd = -rho * 1j * w * np.sum(Phi*n1) # Total drag force
    Cd = 2 * Fd / (rho * U**2 * S)      # Drag coeffient
    
    return f11.real, -f11.imag, Cd
              # /w?     /w**2?

