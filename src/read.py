# -------------- Required packages
import numpy as np
import h5py
import pycstruct
from scipy.interpolate import NearestNDInterpolator

# -------------- Read file functions

def read_hdf5(file):
    # Read hdf5-file:
    h = h5py.File(file, 'r')

    # Internal units:   
    ulength = h['Header'].attrs['UnitLength_in_cm']
    umass   = h['Header'].attrs['UnitMass_in_g']
    uvel    = h['Header'].attrs['UnitVelocity_in_cm_per_s']
    utime   = ulength/uvel
    udens   = umass/(ulength**3)
    uinterg = uvel**2
    iu = {
        'ulength' : ulength,
        'umass'   : umass,
        'uvel'    : uvel,
        'utime'   : utime,
        'udens'   : udens,
        'uinterg' : uinterg
    }

    return h, iu

def read_dump(file, bh=False):
    f = open(file, 'rb')
    
    time = np.fromfile(f, np.float64, 1)
    NSinksAllTasks = np.fromfile(f, np.uint32, 1)
    sinks = {}
        
    fields = ['Pos', 'Vel', 'Accel', 'Mass', 'FormationMass', 'FormationTime', 
              'ID', 'HomeTask', 'Index', 'FormationOrder']
    
    if bh == True:
        fields += ['BlackHole', 'BlackHoleAccretionRadius', 'BlackHoleMass', 'BlackHoleDiskGasMass']
    
    for i in fields:
        sinks[i] = []        
        
    for i in range(NSinksAllTasks[0]):
        struct = pycstruct.StructDef(alignment = 8)
        struct.add('float64', 'Pos', shape=3)
        struct.add('float64', 'Vel', shape=3)
        struct.add('float64', 'Accel', shape=3)
        struct.add('float64', 'Mass')
        struct.add('float64', 'FormationMass')
        struct.add('float64', 'FormationTime')
        struct.add('uint64', 'ID')
        struct.add('uint32', 'HomeTask')
        struct.add('uint32', 'Index')
        struct.add('uint32', 'FormationOrder')
        if bh == True:
            struct.add('uint32', 'BlackHole')
            struct.add('float64', 'BlackHoleAccretionRadius')
            struct.add('float64', 'BlackHoleMass')
            struct.add('float64', 'BlackHoleDiskGasMass')
            
        inbytes = f.read(struct.size())
        data = struct.deserialize(inbytes)
        for field in fields:
            sinks[field] += [data[field]]   
    
    for field in fields:
        sinks[field] = np.array(sinks[field])
    
    f.close()
    
    return time, NSinksAllTasks, sinks


