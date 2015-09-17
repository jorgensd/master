import time; import sys
import warnings
import matplotlib.pyplot as plt, numpy as np
warnings.filterwarnings('ignore')

try:
    a = eval(sys.argv[1])
    b = eval(sys.argv[2])
    N = eval(sys.argv[3])
    geom = sys.argv[4]
    progtype = sys.argv[5]
except IndexError:
    print '\n Usage: %s a b N geom progtype \n' % sys.argv[0]
    print ' a, b     - Major and minor radius for ellipse'
    print '            Height and length for rectangle'
    print ' N        - number of panels'
    print ' geom     - ellipse or rectangle'
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
if geom == 'ellipse':
    from make_grid_ellipse import make_grid_ellipse
    Z = make_grid_ellipse(a, b, N)
else:
    print 'choose ellipse or rectangle'
    sys.exit(0)

tau = np.linspace(0, 1, 50)
m11_list = []; b11_list = []
for t in tau:
    m11, m22, m12, m66, b11 = panel_method(a, b, N, Z, t)
    m11_list.append(m11)
    b11_list.append(b11)
plt.plot(tau, m11_list, tau, b11_list)
plt.legend(['Added mass', 'Damping'])
plt.title('Added mass and damping vs porosity coeff $\\tau$')
plt.xlabel('$\\tau$', fontsize=30)

tau = 0.3; mu = 0.4
A = np.linspace(0, 1, 100)
m11_list = []; b11_list = []
plt.figure()
for A_ in A:
    m11, m22, m12, m66, b11 = panel_method(a, b, N, Z, tau, A_)
    m11_list.append(m11)
    b11_list.append(b11)
KCpor = (1 - tau) * A / (2*a * mu * tau**2)
plt.plot(KCpor, m11_list, KCpor, b11_list)
plt.legend(['Added mass', 'Damping'])
plt.title('Added mass and damping vs KCpor')
plt.xlabel('KCpor', fontsize=17)


plt.show()
