import time; import sys
import warnings
warnings.filterwarnings('ignore')

try:
    a = eval(sys.argv[1])
    b = eval(sys.argv[2])
    N = eval(sys.argv[3])
    progtype = sys.argv[4]
    tau = eval(sys.argv[5])
except IndexError:
    print '\n Usage: %s a b N geom progtype \n' % sys.argv[0]
    print ' a        - Radius'
    print ' N        - number of panels'
    print ' tau      - porosity parameter'
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

start = time.clock()
m11, b11 = panel_method(a, b, N, Z, tau=tau)
end = time.clock()

print '-----------------------------------------------'
print 'Runtime: %.5f s' % (end-start)
print 'Geometry: Cylinder'
print 'Number of panels: %d' % N
print 'Radius: %d' % a
print '-----------------------------------------------'
print '    m11 = ', m11
print '    b11 = ', b11
print '-----------------------------------------------'
