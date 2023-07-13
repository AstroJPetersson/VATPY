# -------------- Required packages
import numpy as np
import h5py
from scipy.interpolate import NearestNDInterpolator

from .read import read_hdf5

# -------------- Frame functions
def number_density(h, iu):
    mp   = 1.6726e-24    #[g]
    
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
    
    numDens['Total'] = np.sum(np.array(list(numDens.values())), axis=0)

    return numDens


def temperature(h, iu):
    mp   = 1.6726e-24    #[g]
    kb   = 1.380658e-16  #[erg K^-1]
    
    if 'ChemicalAbundances' in h['PartType0'].keys():
        # Convert mass densities into number densities:
        numDens = number_density(h, iu)

        # Estimate the temperature of each gas cell:
        interg = h['PartType0']['InternalEnergy'] * iu['uinterg']
        mu     = (numDens['HII'] * 1 + numDens['HI'] * 1 + numDens['H2'] * 2 + numDens['He'] * 4 + numDens['CO'] * 14) / numDens['Total']
        temp   = (2 * mu * mp * interg) / (3 * kb)
    else:
        xHe = 0.1
        mu  = (1 + 4 * xHe)
        interg = h['PartType0']['InternalEnergy'] * iu['uinterg']
        temp   = (2 * mu * mp * interg) / (3 * kb)

    return temp


def density_frame(file, bins=100, axis='z', unit='cgs', cut=None, column=None, box=None):
    # Units:
    pc   = 3.08567758e18                #[cm]
    kpc  = 1e3 * pc
    Myr  = 1e6 * 365.25 * 24 * 60 * 60  #[s]
    
    # Read the data:
    h, iu   = read_hdf5(file=file)
    boxSize = h['Header'].attrs['BoxSize'] * iu['ulength'] / kpc
    time    = h['Header'].attrs['Time'] * iu['utime'] / Myr
    pos     = h['PartType0']['Coordinates'][:] * iu['ulength'] / kpc
    dens    = h['PartType0']['Density'][:] * iu['udens']

    # Generate the grid:
    if box:
        boxMin, boxMax = min(box), max(box) 
    else:
        boxMin, boxMax = 0, boxSize
    X, dX = np.linspace(boxMin, boxMax, bins, retstep=True)
    Y, dY = np.linspace(boxMin, boxMax, bins, retstep=True)
    
    dim = 2
    if cut:
        Z   = cut
        dZ  = 1 / kpc
        dim = 3
    elif column:
        Z, dZ = np.linspace(boxMin, boxMax, bins, retstep=True)
        Z     = Z[(Z > np.min(column)) * (Z < np.max(column))]
    else:
        Z, dZ = np.linspace(boxMin, boxMax, bins, retstep=True)
    XX, YY, ZZ = np.meshgrid(X, Y, Z)

    # Colorbar labels:
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
    
    # Selection of unit on the colorbar:
    if unit != 'cgs':
        numDens = number_density(h, iu)
        dens    = numDens[unit]
        colorbarLabel  = labels[unit]
    else:
        colorbarLabel  = labels['cgs']
   
    # Axis orientation:
    ort = np.empty(np.shape(pos))
    if axis == 'x':
        ort[:,0] = pos[:,1]
        ort[:,1] = pos[:,2]
        ort[:,2] = pos[:,0]
    elif axis == 'y':
        ort[:,0] = pos[:,0]
        ort[:,1] = pos[:,2]
        ort[:,2] = pos[:,1]
    else:
        ort = pos

    # Interpolation:
    interp = NearestNDInterpolator(ort, dens)
    interpDens = interp(XX, YY, ZZ)
    interpDens = np.sum(interpDens * dZ * kpc, axis=2)

    return interpDens, colorbarLabel, time, boxMin, boxMax


def temperature_frame(file, bins=100, axis='z', cut=None, column=None, box=None):
    pc   = 3.08567758e18                #[cm]
    kpc  = 1e3 * pc
    Myr  = 1e6 * 365.25 * 24 * 60 * 60  #[s]

    # Read the data:
    h, iu   = read_hdf5(file=file)
    boxSize = h['Header'].attrs['BoxSize'] * iu['ulength'] / kpc
    time    = h['Header'].attrs['Time'] * iu['utime'] / Myr
    pos     = h['PartType0']['Coordinates'][:] * iu['ulength'] / kpc
    dens    = h['PartType0']['Density'][:] * iu['udens']
    temp    = temperature(h, iu)

    # Coordinates:
    if box:
        boxMin, boxMax = min(box), max(box) 
    else:
        boxMin, boxMax = 0, boxSize
    X, dX = np.linspace(boxMin, boxMax, bins, retstep=True)
    Y, dY = np.linspace(boxMin, boxMax, bins, retstep=True)
    
    dim = 2
    if cut:
        Z = cut
        dZ = 1 / pc
        dim = 3
    elif column:
        Z, dZ = np.linspace(boxMin, boxMax, bins, retstep=True)
        Z = Z[(Z > np.min(column)) * (Z < np.max(column))]
    else:
        Z, dZ = np.linspace(boxMin, boxMax, bins, retstep=True)
    XX, YY, ZZ = np.meshgrid(X, Y, Z)
    
    # Axis orientation:
    ort = np.empty(np.shape(pos))
    if axis == 'x':
        ort[:,0] = pos[:,1]
        ort[:,1] = pos[:,2]
        ort[:,2] = pos[:,0]
    elif axis == 'y':
        ort[:,0] = pos[:,0]
        ort[:,1] = pos[:,2]
        ort[:,2] = pos[:,1]
    else:
        ort = pos
    
    # Interpolation:
    interp_dens = NearestNDInterpolator(ort, dens)
    interp_temp = NearestNDInterpolator(ort, temp)
    interpDens = interp_dens(XX, YY, ZZ)
    interpTemp = interp_temp(XX, YY, ZZ)
    interpTemp = np.sum(interpTemp * interpDens, axis=2) / np.sum(interpDens, axis=2)
    
    return interpTemp, time, boxMin, boxMax


def velocity_frame(file, bins=100, axis='z', cut=None, column=None, box=None):
    pc   = 3.08567758e18                #[cm]
    Myr  = 1e6 * 365.25 * 24 * 60 * 60  #[s]

    # Read the data:
    h, iu   = read_hdf5(file=file)
    boxSize = h['Header'].attrs['BoxSize'] * iu['ulength'] / pc
    time    = h['Header'].attrs['Time'] * iu['utime'] / Myr
    pos     = h['PartType0']['Coordinates'][:] * iu['ulength'] / pc
    vel     = h['PartType0']['Velocities'][:] * iu['uvel'] / 1e5
    dens    = h['PartType0']['Density'][:] * iu['udens']

    # Coordinates:
    if box:
        boxMin, boxMax = min(box), max(box) 
    else:
        boxMin, boxMax = 0, boxSize
    X, dX = np.linspace(boxMin, boxMax, bins, retstep=True)
    Y, dY = np.linspace(boxMin, boxMax, bins, retstep=True)
    
    dim = 2
    if cut:
        Z = cut
        dZ = 1 / pc
        dim = 3
    elif column:
        Z, dZ = np.linspace(boxMin, boxMax, bins, retstep=True)
        Z = Z[(Z > np.min(column)) * (Z < np.max(column))]
    else:
        Z, dZ = np.linspace(boxMin, boxMax, bins, retstep=True)
    XX, YY, ZZ = np.meshgrid(X, Y, Z)
    
    # Axis orientation:
    ort = np.empty(np.shape(pos))
    if axis == 'x':
        ort[:,0] = pos[:,1]
        ort[:,1] = pos[:,2]
        ort[:,2] = pos[:,0]
        vrad     = vel[:,0]
    elif axis == 'y':
        ort[:,0] = pos[:,0]
        ort[:,1] = pos[:,2]
        ort[:,2] = pos[:,1]
        vrad     = vel[:,1]
    else:
        ort  = pos
        vrad = vel[:,2]
    
    # Interpolation:
    interp_dens = NearestNDInterpolator(ort, dens)
    interp_vrad = NearestNDInterpolator(ort, vrad)
    interpDens = interp_dens(XX, YY, ZZ)
    interpVel  = interp_vrad(XX, YY, ZZ)
    interpVel  = np.sum(interpVel * interpDens, axis=2) / np.sum(interpDens, axis=2)

    return interpVel, time, boxMin, boxMax


