import numpy as np
import itertools

def rotate_smallest_volume(pos):
    cov      = np.cov(pos, rowvar=False)
    lam, vec = np.linalg.eigh(cov)
    idx      = np.argsort(lam)[::-1]
    vec      = vec[:, idx]
    pos_out  = np.dot(pos, vec)
    return pos_out

def get_mw(npoly, nbead, mbead, nsolv, msolv):
    return (npoly * nbead * mbead) / (nsolv * msolv + npoly * nbead * mbead)

def get_rho(npoly, nbead, L):
    return npoly * nbead / L ** 3

def gen_data(pos):
    pos   = rotate_smallest_volume(pos) # smallest volume
    vmin  = pos.min(axis=0)
    vmax  = pos.max(axis=0)
    size  = vmax - vmin # e.g., [1.2, 3.1, 2.3] for height, width, and depth 
    
    nbead = len(pos)
    mbead = 5.0         # mass polymer bead
    msolv = 1.0         # mass solvent
    L     = int(5 * nbead ** 0.5 + 0.5) # cubic simulation box side length
    nsolv = 5 * L ** 3  # fixed number density of solvent for neutral buoyancy 
    npoly = L ** 3 / (0.98 / 0.02) / nbead # number of polymer chains 2 wt%
    npoly = int(npoly + 0.5)
    
    lcube = (L / size).astype('int') 
    ncube = np.prod(lcube)  # number of single polymer cubes fit in simulation box
    
    # shrink the original polymer size if necessary
    alpha = 1.0 # shrinking factor
    while ncube < npoly:
        new_size = size / alpha
        lcube    = (L / new_size).astype('int')
        ncube    = np.prod(lcube)
        alpha    = alpha - 0.01
    
    print(f"# Box size: {L}x{L}x{L}")
    print(f"# Shrink r: {alpha}")
    print(f"# Num poly: {npoly}")
    print(f"# Num solv: {nsolv}")
    print(f"# Poly wt : {get_mw(npoly, nbead, mbead, nsolv, msolv):0.6f}")
    print(f"# Poly rho: {get_rho(npoly, nbead, L):0.6f}")
    
    pos_new  = pos - vmin - L / 2
    pos_out  = []
    count = 0
    for i, j, k in itertools.product(range(lcube[0]), range(lcube[1]), range(lcube[2])):
        pos_temp = np.copy(pos_new)
        pos_temp[:, 0] += i * size[0]
        pos_temp[:, 1] += j * size[1]
        pos_temp[:, 2] += k * size[2]
        pos_out.append(pos_temp)
        count += 1
        if count == npoly:
            break
        
    return pos_out, L, nsolv
     
    