'''
Description:

Last updated: 2014-11-12
'''

# -------------- Required packages
import numpy as np
from .read import read_hdf5, read_dump

# -------------- Declare function(s)
def get_black_hole_data(output_dir, n, N, vcr=False):
    # Constants:
    G = 6.67259e-8  # [cm^3 g^-1 s^-2]
    c = 2.998e10  # [cm s^-1]
    mp = 1.6726e-24  # [g]
    kb = 1.380658e-16  # [erg K^-1]
    pc = 3.08567758e18  # [cm]
    kpc = 3.08567758e21  # [cm]
    Myr = 1e6 * 365.25 * 24 * 60 * 60  # [s]
    Msol = 1.9891e33  # [g]
    thom = 6.6525e-25  # [cm^2]

    BlackHoleData = {
        'Time'       : [],
        'MassBH'     : [],
        'MassDisk'   : [],
        'MassReserv' : [],
        'MassSink'   : [],
        'TimeMid'    : [],
        'MdotBH'     : [],
        'MdotEdd'    : [],
        'MdotSink'   : [],
        'AngMom'     : [],
        'CircRadius' : []
    }

    for i in range(n, N+1):
        # Snapshot selection:
        snap = '000'[len(str(i)):] + str(i)
        if i < N:
            print(f'  * Reading data of snapshot {snap}', end='\r')
        else:
            print(f'  * Reading data of snapshot {snap}', end='\n')

        # Read HDF5 & binary dump file:
        h, iu = read_hdf5(f'{output_dir}/snap_{snap}.hdf5')
        dump = read_dump(f'{output_dir}/sink_snap_{snap}', feedback=False,
                         spin=vcr, bh=True, hm=vcr, rcirc=vcr)[2]

        # Append data to the dictionary:
        BlackHoleData['Time'].append(h['Header'].attrs['Time'] 
                                     * iu['utime'] / Myr)
        BlackHoleData['MassBH'].append(dump['BlackHoleMass'][0] 
                                       * iu['umass'] / Msol)
        BlackHoleData['MassDisk'].append(dump['BlackHoleDiskMass'][0] 
                                         * iu['umass'] / Msol)
        BlackHoleData['MassReserv'].append(dump['BlackHoleReservoir'][0] 
                                           * iu['umass'] / Msol)
        BlackHoleData['MassSink'].append(dump['Mass'][0] * iu['umass'] / Msol)
        if vcr == True:
            BlackHoleData['AngMom'].append(dump['AngularMomentum'][0] 
                                           * iu['uangmom'])
            BlackHoleData['CircRadius'].append(dump['BlackHoleCircRadius'][0] 
                                               * iu['ulength'] / pc)

    # Accretion rates:
    dt = (np.array(BlackHoleData['Time'])[1:] 
          - np.array(BlackHoleData['Time'])[:-1])
    mt = ((np.array(BlackHoleData['Time'])[1:] 
           + np.array(BlackHoleData['Time'])[:-1]) / 2)
    BlackHoleData['TimeMid'] = list(mt)
    mdot_sink = ((np.array(BlackHoleData['MassSink'])[1:] 
                  - np.array(BlackHoleData['MassSink'])[:-1]) / (dt * 1e6))
    BlackHoleData['MdotSink'] = list(mdot_sink) 
    mdot_bh = ((np.array(BlackHoleData['MassBH'])[1:] 
                - np.array(BlackHoleData['MassBH'])[:-1]) / (dt * 1e6))
    BlackHoleData['MdotBH'] = list(mdot_bh) 

    # Eddington accretion limit:
    MeanMassBH = ((np.array(BlackHoleData['MassBH'])[:-1] 
                   + np.array(BlackHoleData['MassBH'])[1:]) / 2)
    LimitEdd = (4 * np.pi * G * MeanMassBH * Msol * mp) / (0.1 * thom * c)
    BlackHoleData['MdotEdd'] = list(LimitEdd / (Msol / (365.25 * 24 * 60 * 60)))

    return BlackHoleData

# -------------- End of file
