import time; import sys
import warnings
import matplotlib.pyplot as plt, numpy as np
warnings.filterwarnings('ignore')

try:
    a = eval(sys.argv[1])
    b = eval(sys.argv[2])
    N = eval(sys.argv[3])
    progtype = sys.argv[4]
except IndexError:
    print '\n Usage: %s a b N geom progtype \n' % sys.argv[0]
    print ' a, b     - Major and minor radius for ellipse'
    print ' N        - number of panels'
    print ' progtype - pyvec/pyfort\n'
    sys.exit(0)

if progtype == 'pyvec':
    from panel_method_vectorized import panel_method_vectorized as panel_method
elif progtype == 'pyfort':
    from panel_method import panel_method
else:
    print 'Using vectorized code'
    from panel_method_vectorized import panel_method_vectorized as panel_method

    
# Get grid coordinates
from make_grid_ellipse import make_grid_ellipse
Z = make_grid_ellipse(a, b, N)


tau = np.linspace(0, 1, 50)
m11_list = []; b11_list = []; Cd_list = []
for t in tau:
    m11, b11, Cd = panel_method(a, b, N, Z, t)
    m11_list.append(m11)
    b11_list.append(b11)
    Cd_list.append(Cd)
    
# Plot m11 vs tau
plt.plot(tau, m11_list, tau, b11_list)
plt.legend(['Added mass', 'Damping'])
plt.title('Added mass and damping vs porosity coeff $\\tau$')
plt.xlabel('$\\tau$', fontsize=30)


# Plot Cd vs tau
plt.figure()
plt.plot(tau, Cd_list)
plt.legend(['Drag coefficient'])
plt.title('Drag coefficient vs porosity coeff $\\tau$')
plt.xlabel('$\\tau$', fontsize=30)
plt.ylabel('$C_D$', fontsize=30, rotation=0)


# Plot m11 and b11 vs KCpor
tau = 0.3; mu = 0.4
d = np.linspace(0, 2, 50)
m11_list = []; b11_list = []
plt.figure()
for d_ in d:
    m11, b11, Cd = panel_method(a, b, N, Z, tau=tau, d=d_)
    m11_list.append(m11)
    b11_list.append(b11)
KCpor = (1 - tau) * d / (2*a * mu * tau**2)
plt.plot(KCpor, m11_list, KCpor, b11_list)
plt.legend(['Added mass', 'Damping'], loc='upper left')
plt.title('Added mass and damping vs KCpor')
plt.xlabel('KCpor', fontsize=17)


plt.show()
