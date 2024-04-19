import networkx as nx
import numpy as np


def gen_pos(G, rmin, bdist, bmin, bmax, niter):
    """
    Generate positions for polymer beads in a 3D space based on a given graph structure.

    Parameters:
    -----------
    G : NetworkX graph
        The input graph structure.
    rmin : float
        Minimum separation distance between polymer beads.
    bdist : float
        Ideal bond distance between connected polymer beads.
    bmin : float
        Minimum bond length.
    bmax : float
        Maximum bond length.
    niter : int
        Number of iterations for position generation.

    Returns:
    --------
    numpy.ndarray
        A 2D array containing the positions of polymer beads in a 3D space.

    Notes:
    ------
    This function generates positions for polymer beads in a 3D space based on the given graph structure.
    It uses a combination of force-based methods to compute the positions that satisfy the constraints
    such as separation distance and bond lengths.

    Example:
    --------
    >>> import networkx as nx
    >>> import numpy as np
    >>> G = nx.Graph()
    >>> G.add_edges_from([(0, 1), (1, 2), (2, 3)])
    >>> positions = gen_pos(G, 1.0, 1.0, 0.9, 1.1, 100)
    """
    N    = G.number_of_nodes()
    rcut = rmin / 0.8
    
    # Initial positions using Kamada-Kawai path-length cost-function.
    pos  = nx.kamada_kawai_layout(G, dim=3)
    pos  = np.vstack([pos[key] for key in sorted(pos.keys())])
    adj  = nx.adjacency_matrix(G).toarray()
    
    bonded = adj == 1
    angled = np.matmul(adj, adj) > 0
    
    for step in range(niter):
        print(pos[0])
        r    = pos[:, np.newaxis, :] - pos[np.newaxis, :, :]
        dist = np.sqrt(np.sum(r ** 2, axis=-1))
        
        # Any significant overlap amongst any particles.
        big_overlap = dist < rmin
        big_overlap[bonded] = False
        big_overlap[np.diag_indices(N)] = False
        
        # Any significant stretching and shrinking of the bonds.
        bad_bond = np.logical_or(dist[bonded] < bmin, dist[bonded] > bmax)
        
        if not big_overlap.any() and not bad_bond.any():
            print(f"# No substantial overlaps or bad bonds found after {step} steps...\n\n")
            break
        else:
            if step % 50 == 0:
                print(f"# Significant overlaps or bad bonds detected. Attempting pushoff step {step+1}")
    
        # Force Initialization
        fmag = np.zeros([N, N])
        f    = np.zeros([N, N, 3])
        
        # Compute all pairwise forces using a soft-core potential Uij = A[a+ cos(pi*rij/rc)]
        # Find elements where the distance matrix is less than the cutoff
        
        small_overlap = dist < rcut
        small_overlap[bonded] = False
        small_overlap[np.diag_indices(N)] = False
        
        scale = 0.1
            
        fmag[small_overlap] = scale * np.sin(np.pi * dist[small_overlap] / rcut) / rcut / dist[small_overlap]
        
        # Compute all bonded forces using a harmonic potential Uij = K(d-d0)^2
        fmag[bonded] = - scale * (dist[bonded] - bdist) / dist[bonded]
        
        # TODO: Angle force
        
        # Compute vector force
        f =  r * fmag[:, :, None]
        
        # Get sum of froce on each particle
        fi = np.sum(f, axis=1)
        
        # Evolve positions
        pos += 0.5 * fi
    
    return pos