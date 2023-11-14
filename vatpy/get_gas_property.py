'''
Description: Functions to get various gas properties, e.g. number densities & temperature. 

Last updated: 2023-09-27
'''

import numpy as np


# -------------- Function(s)
def number_density(h, iu):
    # Constants:
    mp   = 1.6726e-24    #[g]
    
    # Gas density:  
    rho = h['PartType0']['Density'] * iu['udens']

    # Chemical abundances:
    X = h['PartType0']['ChemicalAbundances']
    xH2, xHII, xCO = X[:,0], X[:,1], X[:,2]

    # Number densities from chemical abundances:
    xHe = 0.1
    mu = (1 + 4 * xHe)
    nH = rho / (mu * mp)

    numDens = {
        'HII' : xHII * nH, 
        'H2'  : xH2 * nH,
        'HI'  : (1 - xHII - 2*xH2) * nH,
        'CO'  : xCO * nH,
        'He'  : xHe * nH,
        'e'   : xHII * nH
    }
    
    # Total number density:
    numDens['n'] = np.sum(np.array(list(numDens.values())), axis=0)

    return numDens


def temperature(h, iu):
    mp   = 1.6726e-24    #[g]
    kb   = 1.380658e-16  #[erg K^-1]
    
    if 'ChemicalAbundances' in h['PartType0'].keys():
        # Convert mass densities into number densities:
        numDens = number_density(h, iu)

        # Estimate the temperature of each gas cell:
        interg = h['PartType0']['InternalEnergy'] * iu['uinterg']
        mu     = (numDens['HII'] * 1 + numDens['HI'] * 1 + numDens['H2'] * 2 + numDens['He'] * 4 + numDens['CO'] * 14) / numDens['n']
        temp   = (2 * mu * mp * interg) / (3 * kb)
    else:
        xHe = 0.1
        mu  = (1 + 4 * xHe)
        interg = h['PartType0']['InternalEnergy'] * iu['uinterg']
        temp   = (2 * mu * mp * interg) / (3 * kb)

    return temp


# -------------- End of file


