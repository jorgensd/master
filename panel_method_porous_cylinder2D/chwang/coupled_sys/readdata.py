import numpy as np

def readdata(filename):
    """ Reads numerical data and returns two
    arrays containing x and y."""

    infile = open(filename, "r")
    lines = infile.readlines()
    infile.close()
    
    xlist = []
    ylist = []
    for line in lines:
        coor = line.split()
        x = float(coor[1])
        y = float(coor[2])
        xlist.append(x)
        ylist.append(y)
    
    xarr = np.array(xlist)
    yarr = np.array(ylist)
    
    return xarr, yarr
