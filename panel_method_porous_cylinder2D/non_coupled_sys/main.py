import time; import sys
import warnings
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
elif geom == 'rectangle':
    from rectangle_cosspace import rectangle_cosspace
    X, Y = rectangle_cosspace(a, b, N)
else:
    print 'choose ellipse or rectangle'
    sys.exit(0)

start = time.clock()
m11, m22, m12, m66, b11 = panel_method(a, b, N, Z)
end = time.clock()

print '-----------------------------------------------'
print 'Runtime: %.2f s' % (end-start)
print 'Geometry: %s' % geom
print 'Number of panels: %d' % N
if geom == 'ellipse':
    print 'Major radius: %d' % a
    print 'Minor radius: %d' % b
else:
    print 'Height: %d' % a
    print 'Length: %d' % b 
print '-----------------------------------------------'
print '             [%8.5f  %8.5f  %8.5f]' % (m11, abs(m12), 0)
print '    mij =    [%8.5f  %8.5f  %8.5f]' % (abs(m12), m22, 0), ', i,j = 1,2,6.'
print '             [%8.5f  %8.5f  %8.5f]' % (0, 0, m66)
print '-----------------------------------------------'
