from gurobipy import Model, GRB, quicksum, Env
import matplotlib.pyplot as plt
import os

def optimize_node_mis(node_dict, arc_neighborhood, time_limit=600, show_plot=False):
    """
    Optimizes node-based MIS given node dict and arc neighborhood.

    Args:
        node_dict (dict): {node_id (str): [x, y, z]}
        arc_neighborhood (dict): {node_id (str): list of neighboring node_ids}
        time_limit (int): Max optimization time in seconds
        show_plot (bool): Whether to display plot of selected points

    Returns:
        optimal_nodes (list): List of [x, y, z] coordinates selected
        mis_size (int): Number of selected nodes
    """
    options = {
        "WLSACCESSID" : "55c88bbb-cfc0-4d09-b763-e37e0d3bd302",
        "WLSSECRET" : "f21e8941-0ff2-49c2-8500-6883437c1bb2",
        "LICENSEID" : 2677447,
    }
    env = Env(params=options)
    model = Model("node_mis", env=env)
    model.setParam("TimeLimit", time_limit)

    x = model.addVars(node_dict, vtype=GRB.BINARY)

    model.addConstrs(x[i] + (quicksum(x[j] for j in arc_neighborhood[i]) / len(node_dict)) <= 1 for i in node_dict.keys())

    model.setObjective(quicksum(x[i] for i in node_dict.keys()), GRB.MAXIMIZE)
    model.setParam("OutputFlag", 0)
    model.Params.TimeLimit = time_limit
    model.update()

    model.optimize()

    x_vals = model.getAttr("x", x)
    optimal_nodes = [node_dict[i] for i in node_dict if int(x_vals[i]) == 1]

    if show_plot:
        _plot_3d(optimal_nodes)

    return optimal_nodes, len(optimal_nodes), model.MIPGap, model.Runtime

def _plot_3d(points):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for pt in points:
        ax.scatter(pt[0], pt[1], pt[2], s=10, c='green')
    plt.show()
'''
if __name__ == "__main__":
    # Dummy test (replace with real test data)
    dummy_nodes = {
        "0": [1.0, 2.0, 3.0],
        "1": [4.0, 5.0, 6.0],
        "2": [7.0, 8.0, 9.0],
    }
    dummy_arcs = {
        "0": ["1"],
        "1": ["0", "2"],
        "2": ["1"],
    }

    optimal_nodes, size = optimize_node_mis(dummy_nodes, dummy_arcs, show_plot=True)
    print(f"Node MIS size: {size}")
'''