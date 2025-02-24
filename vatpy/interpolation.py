'''
Description: Interpolation of data (e.g. gas physical properties).

Last updated: 2023-09-27
'''

# -------------- Required packages
import numpy as np
from scipy.interpolate import NearestNDInterpolator
from scipy.spatial import KDTree

# -------------- Declare function(s)
def interpolate_to_2d(pos, unit, values, bins, xrange, yrange, zrange, 
                      cut=None, weights=None):
    # Generate a grid:
    X, dX = np.linspace(xrange[0], xrange[1], bins, retstep=True)
    Y, dY = np.linspace(yrange[0], yrange[1], bins, retstep=True)
    mX = (X[:-1] + X[1:]) / 2
    mY = (Y[:-1] + Y[1:]) / 2
    if cut is not None:
        mZ = cut
        dZ = 1 / unit
    else:
        Z, dZ = np.linspace(zrange[0], zrange[1], bins, retstep=True)
        mZ = (Z[:-1] + Z[1:]) / 2
    XX, YY, ZZ = np.meshgrid(mX, mY, mZ)
    coord = np.column_stack((XX.flatten(), YY.flatten(), ZZ.flatten()))

    # Interpolate the data:
    interp = NearestNDInterpolator(pos, values)
    if weights:
        interp_w = NearestNDInterpolator(pos, weights)
    
    # Interpolate onto the grid:
    interpValues = interp(coord)
    if weights:
        interpWeights = interp_w(coord)

    # Sum along the z-axis:
    if np.any(weights):
        H = (np.histogram2d(coord[:,0], coord[:,1], bins=(X, Y), 
                            weights=(interpValues * interpWeights))[0]
             / np.histogram2d(coord[:,0], coord[:,1], bins=(X, Y), 
                              weights=(interpWeights))[0])
    else:
        H = np.histogram2d(coord[:,0], coord[:,1], bins=(X, Y), 
                           weights=(interpValues * dZ * unit))[0]
    H = H.T

    return H


def interpolate_to_2d_kdtree(pos, unit, values, bins, xrange, yrange, zrange, 
                             cut=None, weights=None):
    # Tree construction:
    tree = KDTree(pos, leafsize=10)

    # Generate a grid:
    X, dX = np.linspace(xrange[0], xrange[1], bins, retstep=True)
    Y, dY = np.linspace(yrange[0], yrange[1], bins, retstep=True)
    mX = (X[:-1] + X[1:]) / 2
    mY = (Y[:-1] + Y[1:]) / 2
    if cut is not None:
        mZ = cut
        dZ = 1 / unit
    else:
        Z, dZ = np.linspace(zrange[0], zrange[1], bins, retstep=True)
        mZ = (Z[:-1] + Z[1:]) / 2
    XX, YY, ZZ = np.meshgrid(mX, mY, mZ)
    coord = np.column_stack((XX.flatten(), YY.flatten(), ZZ.flatten()))

    # Find Indices:
    q = tree.query(coord)
    
    # Sum along the z-axis:
    if np.any(weights):
        H = (np.histogram2d(coord[:,0], coord[:,1], bins=(X, Y), 
                            weights=(values[q[1]] * weights[q[1]]))[0]
             / np.histogram2d(coord[:,0], coord[:,1], bins=(X, Y), 
                              weights=(weights[q[1]]))[0])
    else:
        H = np.histogram2d(coord[:,0], coord[:,1], bins=(X, Y), 
                           weights=(values[q[1]] * dZ * unit))[0]
    H = H.T

    return H

# -------------- End of file

