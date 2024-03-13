import numpy as np
from graph2lammps.convert import gen_pos

def write_pdb(G, pos, fname):
    fid = open(fname, "w")
    fid.write("COMPND    UNNAMED\n")
    fid.write("AUTHOR    GENEATED BY graph2lammps.py\n")
  
    # WRITE HETATM ENTRIES.
    for i, (node, data) in enumerate(G.nodes(data=True)):
        fid.write("{:<6s}{:>5d}  {:<3s}{:1s}{:>3s}  {:>4d}{:<1s}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}      {:>4s} {:<s}\n".\
                format("ATOM", data["type"], data["name"], " ", "UNL", i+1, " ", \
                pos[i, 0], pos[i, 1], pos[i, 2], 1.00, 0.00, "", data["name"]))
        
    # WRITE CONECT ENTRIES.
    for edge in G.edges():
        source, target = edge
        fid.write("CONECT {:>5d} {:>5d}\n".format(source+1, target+1))
        
        
def write_lammps(G, pos, fname):
    nnodes = G.number_of_nodes()
    nedges = G.number_of_edges()
    
    fid = open(fname, "w")
    fid.write('LAMMPS GENEATED BY graph2lammps.py\n\n')
    fid.write(f'{nnodes:<6} atoms\n')
    fid.write(f'{nedges:<6} bonds\n')
    
    # Node/bead/atom types.
    node_types  = np.array([G.nodes[node].get('type', 1) for node in G.nodes()])
    node_types  = np.unique(node_types)
    nnode_types = len(node_types)
    fid.write(f'{nnode_types} atom types\n')
    
    # Edge/bond types.
    edge_types  = np.array([G.edges[edge].get('type', 1) for edge in G.edges()])
    edge_types  = np.unique(edge_types)
    nedge_types = len(edge_types)
    fid.write(f'{nedge_types} bond types\n\n')
    
    # Simulation box.
    # TODO: User-input box config; transform pos.
    min_x, min_y, min_z = np.min(pos, axis=0)
    max_x, max_y, max_z = np.max(pos, axis=0)
    fid.write(f'{min_x-5:<7.4f} {max_x+5:<7.4f} xlo xhi\n')
    fid.write(f'{min_y-5:<7.4f} {max_y+5:<7.4f} ylo yhi\n')
    fid.write(f'{min_z-5:<7.4f} {max_z+5:<7.4f} zlo zhi\n\n')
    
    # Bead masses.
    fid.write('Masses\n\n')
    for key, num in G.graph['mass'].items():
        fid.write(f'{key:<3} {num:<}\n')
    fid.write('\n')
    
    # Bead positions.
    # TODO: a more general case.
    fid.write('Atoms\n\n')
    for i in range(len(pos)):
        n_type = G.nodes[i]['type']
        fid.write(f'{i+1:<7} 1 {n_type} 0 {pos[i, 0]:<7.4f} {pos[i, 1]:<7.4f} {pos[i, 2]:<7.4f}\n')
        
    # Bond info.
    fid.write('\n\nBonds\n\n')
    for i, edge in enumerate(G.edges()):
        source, target = edge
        e_type = G.edges[edge].get('type', 1)
        fid.write(f'{i+1:<4} {e_type:<4} {source+1:<4} {target+1:>4}\n')
        
    # TODO: angle, dihedral