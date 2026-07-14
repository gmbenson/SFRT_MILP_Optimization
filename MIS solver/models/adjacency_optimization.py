from gurobipy import Model, GRB, quicksum, Env
import numpy as np
import matplotlib.pyplot as plt

def optimize_adjacency_mis(node_dict, adjacency_matrix, time_limit=600, show_plot=False):
    """
    Runs MIS optimization using an adjacency matrix.

    Args:
        node_dict (dict): {node_id (str): [x, y, z]}
        adjacency_matrix (2D list or np.array): adjacency[i][j] = 1 if nodes are connected
        time_limit (int): Gurobi solver time limit in seconds
        show_plot (bool): Whether to show 3D plot

    Returns:
        optimal_nodes (list): Selected nodes as [x, y, z] coordinates
        mis_size (int): Size of the MIS
    """
    options = {
        "WLSACCESSID" : "55c88bbb-cfc0-4d09-b763-e37e0d3bd302",
        "WLSSECRET" : "f21e8941-0ff2-49c2-8500-6883437c1bb2",
        "LICENSEID" : 2677447,

    }
    env = Env(params=options)
    model = Model("adjacency_mis", env=env)
    model.setParam("OutputFlag", 0)
    model.setParam("TimeLimit", time_limit)

    nodes = list(node_dict.keys())  # Assume string keys
    x = model.addVars(nodes, vtype=GRB.BINARY)

    a = np.identity(len(nodes), dtype=int) - np.array(adjacency_matrix)

    model.setObjective(
        quicksum(x[i] * a[int(i)][int(j)] * x[j] for i in nodes for j in nodes),
        GRB.MAXIMIZE
    )

    model.optimize()
    x_vals = model.getAttr("x", x)

    optimal_nodes = [node_dict[i] for i in nodes if int(x_vals[i]) == 1]

    if show_plot:
        _plot_3d(optimal_nodes)

    return optimal_nodes, len(optimal_nodes), model.MIPGap, model.Runtime

def _plot_3d(points):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for pt in points:
        ax.scatter(pt[0], pt[1], pt[2], s=10, c='blue')
    plt.show()


# Optional CLI for debugging
if __name__ == "__main__":
    import sys

    node_file = "../data/node_data.pkl"
    adjacency_file = "../data/adjacency_matrix.pkl"

    optimal_nodes, size = optimize_adjacency_mis(node_file, adjacency_file, show_plot=True)
    print(f"Size of adjacency MIS: {size}")
