# -------------- Required packages
import numpy as np
import h5py
from scipy.interpolate import NearestNDInterpolator

from .read import read_hdf5

# -------------- Frame functions

def number_density(h, iu):
    # Gas density:  
    rho = h['PartType0']['Density'] * iu['udens']

    # Chemical abundances:
    X = h['PartType0']['ChemicalAbundances']
    xH2, xHII, xCO = X[:,0], X[:,1], X[:,2]

    # Number density of hydrogen:
    xHe = 0.1
    mu = (1 + 4 * xHe)
    nH = rho / (mu * mp)

    # More number densities:
    numDens = {
        'HII' : xHII * nH, 
        'H2'  : xH2 * nH,
        'HI'  : (1 - xHII - 2*xH2) * nH,
        'CO'  : xCO * nH,
        'He'  : xHe * nH,
        'e'   : xHII * nH
    }
    
    numDens['Total'] = np.sum(np.array(list(n.values())), axis=0)

    return numDens


def density_frame(file, bins=100, lookDownAxis='z', unit='cgs', zSlice=None, zColumn=None):
    pc   = 3.08567758e18                #[cm]
    Myr  = 1e6 * 365.25 * 24 * 60 * 60  #[s]
    
    # Read data:
    h, iu   = read_hdf5(file=file)
    boxSize = h['Header'].attrs['BoxSize'] * iu['ulength'] / pc
    time    = h['Header'].attrs['Time'] * iu['utime'] / Myr
    pos     = h['PartType0']['Coordinates'][:] * iu['ulength'] / pc
    dens    = h['PartType0']['Density'][:] * iu['udens']
    
    # Coordinates:
    X, dX = np.linspace(0, boxSize, bins, retstep=True)
    Y, dY = np.linspace(0, boxSize, bins, retstep=True)
    dim = 2
    if zSlice != None:
        Z = zSlice
        dZ = 1 / pc
        dim = 3
    elif zColumn != None:
        Z, dZ = np.linspace(0, boxSize, bins, retstep=True)
        Z = Z[(Z > np.min(zColumn)) * (Z < np.max(zColumn))]
    else:
        Z, dZ = np.linspace(0, boxSize, bins, retstep=True)
    XX, YY, ZZ = np.meshgrid(X, Y, Z)
    
    # Surface density labels:
    labels = {
        'HII' : r'$\log_{10}(\Sigma_{\mathrm{HII}} \ [\mathrm{cm}^{-%d}])$' % dim,
        'H2'  : r'$\log_{10}(\Sigma_{\mathrm{H}_2} \ [\mathrm{cm}^{-%d}])$' % dim,
        'HI'  : r'$\log_{10}(\Sigma_{\mathrm{HI}} \ [\mathrm{cm}^{-%d}])$' % dim,
        'CO'  : r'$\log_{10}(\Sigma_{\mathrm{CO}} \ [\mathrm{cm}^{-%d}])$' % dim,  
        'He'  : r'$\log_{10}(\Sigma_{\mathrm{He}} \ [\mathrm{cm}^{-%d}])$'% dim,  
        'e'   : r'$\log_{10}(\Sigma_{e^{-}} \ [\mathrm{cm}^{-%d}])$' % dim,
        'n'   : r'$\log_{10}(\Sigma_\mathrm{Gas} \ [\mathrm{cm}^{-%d}])$' % dim,
        'cgs' : r'$\log_{10}(\Sigma_\mathrm{Gas} \ [\mathrm{g} \ \mathrm{cm}^{-%d}])$' % dim
    }
    
    # Units:
    if unit != 'cgs':
        numDens = number_density(h, iu)
        dens    = numDens[unit]
        colorbarLabel  = labels[unit]
    else:
        colorbarLabel  = labels['cgs']
    
    # Interpolation:
    axis = {'x' : 0, 'y' : 1, 'z' : 2}
    interp = NearestNDInterpolator(pos, dens)
    interpDens = interp(XX, YY, ZZ)
    interpDens = np.sum(interpDens * dZ * pc, axis=axis[lookDownAxis])

    return interpDens, colorbarLabel, time, boxSize


