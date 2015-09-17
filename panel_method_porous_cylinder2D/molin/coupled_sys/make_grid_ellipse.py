import numpy as np

def make_grid_ellipse(r_a, r_b, N):
    # r_a - major radius
    # r_b - minor radius
    # N - number of panels
    #
    # returns the the complex coordinates

    x = r_a*np.cos( np.linspace(0, 2*np.pi, N+1) )
    y = r_b*np.sin( np.linspace(0, 2*np.pi, N+1) )
    z = x + 1j*y
    
    return z
