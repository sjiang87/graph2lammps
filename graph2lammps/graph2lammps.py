#!/usr/bin/env python
__author__ = "Shengli Jiang"
__credits__ = ["Shengli Jiang", "Michael A. Webb"]
__email__ = "sj0161@princeton.edu"

"""
This program (and other distributed files) enable the facile
generation of LAMMPS polymer data based on networkx graphs.
The position function is adapted from Prof. Webb's polymerize.py.
"""

from graph2lammps.convert import gen_pos
from graph2lammps.convert import write_lammps


def graph2lammps(G, rmin=0.5, bdist=1.0, bmin=0.9, bmax=1.1, niter=1000, fname='./sys.data'):
    """
    Convert polymer graph structure to LAMMPS input file.

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
    None
    """

    print("# Wring out polymer structure to ",fname," ...")
    
    pos = gen_pos(G, 
                  rmin=rmin, 
                  bdist=bdist, 
                  bmin=bmin, 
                  bmax=bmax, 
                  niter=niter)
    
    write_lammps(G, 
                 pos, 
                 fname=fname)
    
    return None
    

