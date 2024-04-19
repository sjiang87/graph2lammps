import numpy as np


def read_pdb(file):
    pos = []
    with open(file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                tokens = line.split()
                x, y, z = float(tokens[6]), float(tokens[7]), float(tokens[8])
                pos.append([x, y, z])
            elif line.startswith('CONECT'):
                break  # Assuming atom positions come before connectivity information
    return np.array(pos)