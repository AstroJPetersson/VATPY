'''
Description: Functions to get images for the VATPY CLI. 

Last updated: 2024-07-10
'''

# -------------- Required packages
import numpy as np
from scipy.interpolate import NearestNDInterpolator
from scipy.spatial.transform import Rotation
from .read import read_hdf5
from .get_gas_property import temperature

# -------------- Declare function(s)
def get_gas_density_image_cli(h, iu, axis='z', rotate=0, bins=100, xrange=None, yrange=None, 
                              zrange=None, quantity='gas', ulength='kpc'):
    # Unit selection:
    if ulength == 'kpc':
        ulength = 3.08567758e21
    else:
        ulength = 1

    # Get the necessary data: 
    time    = h['Header'].attrs['Time'] * iu['utime']
    boxsize = h['Header'].attrs['BoxSize'] * iu['ulength'] / ulength
    pos     = h['PartType0']['Coordinates'] * iu['ulength'] / ulength
    dens    = h['PartType0']['Density'] * iu['udens']

    # Quantity selection:
    if (quantity != 'gas'):
        num  = number_density(h, iu)
        dens = numDens[quantity]

    # Rotation:
    if rotate != 0:
        rotation = Rotation.from_euler(axis, rotate, degrees=True)
        pos = rotation.apply(pos - boxsize/2)
        pos += boxsize/2
    
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

def get_gas_temperature_image_cli(h, iu, axis='z', rotate=0, bins=100, xrange=None, yrange=None, 
                                  zrange=None, ulength='kpc'):
    # Unit selection:
    if ulength == 'kpc':
        ulength = 3.08567758e21
    else:
        ulength = 1
    
    # Get the necessary data: 
    time    = h['Header'].attrs['Time'] * iu['utime']
    boxsize = h['Header'].attrs['BoxSize'] * iu['ulength'] / ulength
    pos     = h['PartType0']['Coordinates'] * iu['ulength'] / ulength
    dens    = h['PartType0']['Density'] * iu['udens']
    temp    = temperature(h, iu)

    # Rotation:
    if rotate != 0:
        rotation = Rotation.from_euler(axis, rotate, degrees=True)
        pos = rotation.apply(pos - boxsize/2)
        pos += boxsize/2
    
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


