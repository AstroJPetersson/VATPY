'''
Description: Interpolation of data (e.g. gas physical properties) onto a grid.

Last updated: 2023-09-27
'''

# -------------- Required packages
import numpy as np
from scipy.interpolate import NearestNDInterpolator

# -------------- Function(s)
def interpolate_to_2d(pos, unit, values, bins, xrange, yrange, zrange, cut=None):
    # Generate the grid:
    X, dX = np.linspace(xrange[0], xrange[1], bins, retstep=True)
    Y, dY = np.linspace(xrange[0], yrange[1], bins, retstep=True)
    if cut:
        Z   = cut
        dZ  = 1
    else:
        Z, dZ = np.linspace(zrange[0], zrange[1], bins, retstep=True)
    XX, YY, ZZ = np.meshgrid(X, Y, Z)
    coord = np.column_stack((XX.flatten(), YY.flatten(), ZZ.flatten()))

    # Interpolate the data:
    interp = NearestNDInterpolator(pos, values)
    interpValues = interp(coord)

    # Now make a 2d histogram:
    H = np.histogram2d(coord[:,0], coord[:,1], bins=(X, Y), weights=(interpValues * dZ * unit))[0]
    H = H.T

    return H

# -------------- End of file

