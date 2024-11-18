import numpy as np
import itertools


def rotate_smallest_volume(pos):
    cov = np.cov(pos, rowvar=False)
    lam, vec = np.linalg.eigh(cov)
    idx = np.argsort(lam)[::-1]
    vec = vec[:, idx]
    pos_out = np.dot(pos, vec)
    return pos_out


def uniform_sample(vmin, vmax, n):
    ranges = vmax - vmin
    total_points = np.prod(ranges)
    n = min(n, total_points)

    # Generate uniform samples along each axis
    samples = [
        np.linspace(vmin[i], vmax[i], int(ranges[i]), endpoint=False)
        for i in range(len(vmin))
    ]

    indices = np.random.choice(total_points, n, replace=False)

    points = np.array(
        [samples[i][indices % int(ranges[i])] for i in range(len(vmin))]).T

    return points.astype('int')


def get_mw(npoly, nbead, mbead, nsolv, msolv):
    return (npoly * nbead * mbead) / (nsolv * msolv + npoly * nbead * mbead)


def get_rho(npoly, nbead, L):
    return npoly * nbead / L**3


def gen_data(pos, args):
    pos = rotate_smallest_volume(pos)  # smallest volume

    vmin = pos.min(axis=0)
    vmax = pos.max(axis=0)
    size = vmax - vmin + .2  # e.g., [1.2, 3.1, 2.3] for height, width, and depth, add 0.5, so they don't touch

    nbead = len(pos)
    mbead = 5.0  # mass polymer bead
    msolv = 1.0  # mass solvent
    L = 75  # int(5 * nbead ** 0.5 + 0.5) # cubic simulation box side length
    nsolv = 5 * L**3  # fixed number density of solvent for neutral buoyancy
    npoly = L**3 / (
        (1 - args.rho) / args.rho) / nbead  # number of polymer chains 2 wt%
    npoly = int(npoly + 0.5)

    lcube = (L / size).astype('int')
    ncube = np.prod(
        lcube)  # number of single polymer cubes fit in simulation box
    # shrink the original polymer size if necessary
    alpha = 1.0  # shrinking factor
    new_size = np.copy(size)
    while ncube < npoly and alpha > 0.02:
        new_size = size * alpha
        lcube = (L / new_size).astype('int')
        ncube = np.prod(lcube)
        alpha = alpha - 0.01

    if alpha == 1.0:
        new_size = size
    else:
        alpha = alpha + 0.01

    if alpha < 0.02:
        raise RuntimeError("# Shrinking factor too small ...")

    print(f"# Box size: {L}x{L}x{L}")
    print(f"# Cube    : {lcube[0]}x{lcube[1]}x{lcube[2]}")
    print(f"# N Cube  : {ncube}")
    print(f"# Shrink r: {alpha}")
    print(f"# Num poly: {npoly}")
    print(f"# Num solv: {nsolv}")
    print(f"# Poly wt : {get_mw(npoly, nbead, mbead, nsolv, msolv):0.6f}")
    print(f"# Poly rho: {get_rho(npoly, nbead, L):0.6f}")

    pos_new = np.copy(pos)
    pos_new = pos * alpha
    vmin_new = pos_new.min(axis=0)
    pos_new = pos_new - vmin_new
    pos_new = pos_new - L / 2
    pos_out = []
    count = 0

    for i, j, k in itertools.product(range(lcube[0]), range(lcube[1]),
                                     range(lcube[2])):
        pos_temp = np.copy(pos_new)
        pos_temp[:, 0] += i * new_size[0]
        pos_temp[:, 1] += j * new_size[1]
        pos_temp[:, 2] += k * new_size[2]
        pos_temp += np.random.random((pos.shape)) * 0.001
        pos_out.append(pos_temp)
        count += 1
        if count == npoly:
            break

    # points = uniform_sample(vmin=np.array([0, 0, 0]), vmax=lcube, n=npoly)
    # for p in points:
    #     pos_temp = np.copy(pos_new)
    #     pos_temp[:, 0] += p[0] * size[0]
    #     pos_temp[:, 1] += p[1] * size[1]
    #     pos_temp[:, 2] += p[2] * size[2]
    #     pos_out.append(pos_temp)

    return pos_out, L, nsolv, npoly
