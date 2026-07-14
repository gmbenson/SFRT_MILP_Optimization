import math
import json
import pickle

def conversion(z, point, nodes):
    index = 0
    i = 0
    while i <= z:  # and i <= len(nodes) - 1:
        j = 0
        while j <= point and j <= len(nodes[i]) - 1:
            index += 1
            j += 1
        i += 1
    return index

def create_edges(nodes, nodes_dict, min_dist=30, z_space=2.5):  # z_space must be the space between z planes
    edges = []
    z_list = []

    for i in range(len(nodes)):
        if len(nodes[i]) != 0:
            z_list.append(nodes[i][0][2])

    for z in range(len(z_list)):
        z_min = z_list[z] - min_dist
        z_max = z_list[z] + min_dist
        z_min_holder = z_list[0]
        z_max_holder = z_list[-1]

        for i in z_list:
            if i <= z_min and abs(i - z_min) < abs(z_min_holder - z_min):
                z_min_holder = i
            if i >= z_max and i - z_max < z_max_holder - z_max:
                z_max_holder = i

        z_vals = []
        while z_min_holder <= z_max_holder:
            z_vals.append(z_min_holder)
            z_min_holder += z_space

        for i in range(len(nodes[z])):
            for j in range(len(z_vals)):
                for k in range(len(nodes[j])):
                    if math.dist(nodes[z][i], nodes[j][k]) <= min_dist and nodes[z][i] != nodes[j][k]:
                        edges.append([conversion(z, i, nodes), conversion(j, k, nodes)])
        print(f"{z} out of {len(z_list)} slices done")
    return edges

def create_arcs(nodes, min_dist=30):
    arcs = []
    node_arcs = {}
    matrix = []
    neighborhood = {key: [] for key in nodes}
    num_arcs = 0
    for i in nodes.keys():
        holder = []
        matrix_holder = []
        for j in nodes.keys():
            if math.dist(nodes[i], nodes[j]) <= min_dist and i != j:
                arcs.append([i, j])
                holder.append(j)
                matrix_holder.append(1)
                num_arcs += 1
            elif min_dist < math.dist(nodes[i], nodes[j]) <= min_dist + 10 and i != j:
                neighborhood[i].append(j)
                matrix_holder.append(0)
            else:
                matrix_holder.append(0)
        node_arcs[i] = holder
        matrix.append(matrix_holder)
    return arcs, node_arcs, matrix, neighborhood, num_arcs

def process_nodes(node_list, output_prefix, min_dist=30, z_space=2.5):
    """
    Main function to process node list, create arcs, and save data.
    Returns (num_nodes, num_arcs)
    """
    node_dict = {}
    index = 0
    for point in node_list:
        node_dict[str(index)] = point
        index += 1

    num_nodes = len(node_dict)
    print(f"Loaded {num_nodes} nodes")

    arcs, node_arcs, matrix, neighborhood, num_arcs = create_arcs(node_dict, min_dist=min_dist)
    print(f"Created {num_arcs} arcs")

    # Save outputs with output_prefix
    with open(f'{output_prefix}_node_data.pkl', 'wb') as outfile:
        pickle.dump(node_dict, outfile)

    with open(f'{output_prefix}_arc_data.pkl', 'wb') as outfile:
        pickle.dump(arcs, outfile)

    with open(f'{output_prefix}_adjacency_matrix.pkl', 'wb') as outfile:
        pickle.dump(matrix, outfile)

    with open(f'{output_prefix}_arc_list.pkl', 'wb') as outfile:
        pickle.dump(node_arcs, outfile)

    with open(f'{output_prefix}_neighborhood.pkl', 'wb') as outfile:
        pickle.dump(neighborhood, outfile)

    # Inside helpers/parameter_writer.py
    return {
        "nodes": node_dict,
        "arcs": arcs,
        "adjacency": matrix,
        "num_nodes": len(node_dict),
        "num_arcs": num_arcs,
        "neighborhood": node_arcs
    }


if __name__ == '__main__':
    # Example usage when running parameter_writer.py standalone:
    with open("../data/candidate_points008.json", 'r') as openfile:
        candidate_points = json.load(openfile)

    num_nodes, num_arcs = process_nodes(candidate_points, output_prefix='../data/param_output008')
    print(f"Graph size: {num_nodes} nodes, {num_arcs} edges")
