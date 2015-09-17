import numpy as np

def rectangle_cosspace(a, b, N):
    # make_grid generates a grid for use with the
    # panel method, using cosine spacing.
    #
    # a - rectangle height
    # b - rectangle length
    # N - number of panels
    #
    # returns the x- and y-coordinates

    x = np.zeros(N+1)
    y = np.zeros(N+1)
    a = a/2. # distance from the origin
    b = b/2.

    # This method requires N to be a multiplum of 4
    N = N/4 * 4
    beta = np.linspace(0, np.pi, N/4. + 2)
    # Exlude the endpoints to avoid sharp edges
    alfa = beta[1:-1]
    
    # Right side vertical, x = a
    y_rs = -a * np.cos(alfa)
    x_rs = b * np.ones(len(y_rs))

    # Upper side horisontal, y = a
    x_us = b * np.cos(alfa)
    y_us = a * np.ones(len(x_us)) 

    # Left side vertical, x = -b
    y_ls = a * np.cos(alfa)
    x_ls = -b * np.ones(len(y_ls))

    # Bottom side horisontal, y = -a
    x_bs = -b * np.cos(alfa)
    y_bs = -a * np.ones(len(x_bs))
    
    for x_tmp in [[x_us,y_us],[x_ls,y_ls],[x_bs,y_bs],[x_rs[0],y_rs[0]]]:
        x_rs = np.append(x_rs, x_tmp[0])
        y_rs = np.append(y_rs, x_tmp[1])

    return x_rs, y_rs

if __name__ == '__main__':
    import sys, matplotlib.pyplot as plt
    a = eval(sys.argv[1]); b = eval(sys.argv[2])
    N = int(sys.argv[3])
    x, y = rectangle_cosspace(a, b, N)
    plt.plot(x,y, 'o-')
    plt.axis([-a*0.7, a*0.7, -b*0.7, b*0.7])
    plt.show()
