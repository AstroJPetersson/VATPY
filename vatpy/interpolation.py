'''
Description: Interpolation of e.g. gas properties onto a grid.

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

    # Interpolate the data:
    interp = NearestNDInterpolator(pos, values)
    interpValues = interp(XX, YY, ZZ)
    interpValues = np.sum(interpValues * dZ * unit, axis=2)

    return interpValues


# -------------- End of file


