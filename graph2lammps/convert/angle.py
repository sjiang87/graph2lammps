def angle_index(G):
    indices = []
    for node in G.nodes():
        neighbors = set(G.neighbors(node))
        for neighbor1 in neighbors:
            for neighbor2 in neighbors - {neighbor1}:
                indices.append([node, neighbor1, neighbor2])
    return indices