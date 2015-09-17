import numpy as np, sys
import cmath as cm
from math import log, sqrt

def panel_method_vectorized(a, b, N, Z):
    
    cz = 0.5 * (Z[1:] + Z[:-1])           # complex midpoint
    cx = cz.real                          # x-coor panel midpoint
    cy = cz.imag                          # y-coor panel midpoint
    L = abs(Z[1:] - Z[:-1])               # Panel length
    n1 = (-Z.imag[1:] + Z.imag[:-1]) / L  # x-comp normal vector
    n2 = ( Z.real[1:] - Z.real[:-1]) / L  # y-comp normal vector
    x = Z.real
    y = Z.imag
    
    # Matrix coefficients P and q
    P = np.zeros([N, N])
    q = np.zeros([N, N])

    #
    def TpG(x0, x1, xm, y0, y1, ym, t):
        r = np.sqrt((0.5*(x0*(1-t)+x1*(1+t))-xm)**2+(0.5*(y0*(1-t)+y1*(1+t))-ym)**2)
        v = 0.5 * np.sqrt((x1 - x0)**2 + (y1 - y0)**2)
        return np.log(r)*v

    #
    def TpG2(x0, x1, xm, y0, y1, ym, t):
        r = np.sqrt((0.5*(x0*(1-t)+x1*(1+t))-xm)**2+(0.5*(y0*(1-t)+y1*(1+t))-ym)**2)
        v = 0.5 * np.sqrt((x1 - x0)**2 + (y1 - y0)**2)
        return v / r
    
    # Weights and points for gauss quad
    from readdata import readdata
    w, p = readdata('weightAndPoints.txt')
        
    # Opening angle given by the law of cosines
    for i in xrange(N):
        b = abs(cz[i] - Z[1:])
        c = abs(cz[i] - Z[:-1])
        #P[i] = -np.arccos((b**2 + c**2 - L**2) / (2*b*c))
        for (wi, pi) in zip(w, p):
                q[i] += wi * TpG(x[:-1], x[1:], cx[i], y[:-1], y[1:], cy[i], pi)
                P[i] -= wi * TpG2(x[:-1], x[1:], cx[i], y[:-1], y[1:], cy[i], pi)

    P[np.isnan(P)] = 0 # Needed for rectangle
    #np.fill_diagonal(P, -np.pi) # phi = -pi, for i = j
    
    # RHS
    Q = np.transpose([np.dot(q,n1),np.dot(q,n2),np.dot(q,(cx*n2 - cy*n1))])

    # Velocity potential for each panel
    phi_i = np.linalg.solve(P, Q)
    
    # Added mass
    m11 = np.sum( phi_i[:, 0] * n1 * L )       # m11 for the whole body
    m22 = np.sum( phi_i[:, 1] * n2 * L )
    m12 = np.sum( phi_i[:, 0] * n2 * L )       # Cross coupling
    m66 = np.sum( phi_i[:, 2] * (cx*n2 - cy*n1) * L )

    return m11, m22, m12, m66
