'''
Description: Functions to get images for the VATPY GUI. 

Last updated: 2023-09-27
'''

# -------------- Required packages
import numpy as np
from scipy.interpolate import NearestNDInterpolator
from .read import read_hdf5
from .get_gas_property import temperature

# -------------- Declare function(s)
def get_gas_density_image(file, bins, xrange, yrange, zrange, ulength):
    # Read the data: 
    h, iu   = read_hdf5(file)
    time    = h['Header'].attrs['Time'] * iu['utime']
    boxsize = h['Header'].attrs['BoxSize'] * iu['ulength'] / ulength
    pos     = h['PartType0']['Coordinates'] * iu['ulength'] / ulength
    dens    = h['PartType0']['Density'] * iu['udens']

    # Coordinate ranges:
    if not (xrange):
        xrange = (0, boxsize)
    if not (yrange):
        yrange = (0, boxsize)
    if not (zrange):
        zrange = (0, boxsize)

    # Generate a grid:
    X, dX = np.linspace(xrange[0], xrange[1], bins, retstep=True)
    Y, dY = np.linspace(yrange[0], yrange[1], bins, retstep=True)
    Z, dZ = np.linspace(zrange[0], zrange[1], bins, retstep=True)
    XX, YY, ZZ = np.meshgrid(X, Y, Z)
    coord = np.column_stack((XX.flatten(), YY.flatten(), ZZ.flatten()))

    # Interpolate the data:
    interpolator = NearestNDInterpolator(pos, dens)

    # Interpolate onto the grid:
    interp = interpolator(coord)
    H = np.histogram2d(coord[:,0], coord[:,1], bins=bins, weights=(interp * dZ * ulength))[0]
    H = np.log10(H.T)

    return H, X, Y

def get_gas_temperature_image(file, bins, xrange, yrange, zrange, ulength):
    # Read the data: 
    h, iu   = read_hdf5(file)
    time    = h['Header'].attrs['Time'] * iu['utime']
    boxsize = h['Header'].attrs['BoxSize'] * iu['ulength'] / ulength
    pos     = h['PartType0']['Coordinates'] * iu['ulength'] / ulength
    dens    = h['PartType0']['Density'] * iu['udens']
    temp    = temperature(h, iu)

    # Coordinate ranges:
    if not (xrange):
        xrange = (0, boxsize)
    if not (yrange):
        yrange = (0, boxsize)
    if not (zrange):
        zrange = (0, boxsize)

    # Generate a grid:
    X, dX = np.linspace(xrange[0], xrange[1], bins, retstep=True)
    Y, dY = np.linspace(yrange[0], yrange[1], bins, retstep=True)
    Z, dZ = np.linspace(zrange[0], zrange[1], bins, retstep=True)
    XX, YY, ZZ = np.meshgrid(X, Y, Z)
    coord = np.column_stack((XX.flatten(), YY.flatten(), ZZ.flatten()))

    # Interpolate the data:
    interpolator_dens = NearestNDInterpolator(pos, dens)
    interpolator_temp = NearestNDInterpolator(pos, temp)

    # Interpolate onto the grid:
    interp_dens = interpolator_dens(coord)
    interp_temp = interpolator_temp(coord)
    H = (np.histogram2d(coord[:,0], coord[:,1], bins=bins, weights=(interp_temp * interp_dens))[0]
         / np.histogram2d(coord[:,0], coord[:,1], bins=bins, weights=(interp_dens))[0])
    H = np.log10(H.T)

    return H, X, Y

# -------------- End of file


