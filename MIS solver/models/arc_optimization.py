import math
import pickle
import csv
from gurobipy import *
import matplotlib.pyplot as plt
import plotly.express as px


def optimize_arc_mis(node_dict, arcs, time_limit=600, show_plot=False):
    """
    Solves Maximum Independent Set problem using arc constraints (pairwise non-adjacency).

    Args:
        node_dict (dict): {str: [x, y, z]} node ID → coordinate
        arcs (list): list of [i, j] edges where i, j are str IDs
        time_limit (int): max time in seconds
        show_plot (bool): optionally show plots

    Returns:
        optimal_nodes (list): list of [x, y, z] points in the MIS
        mis_size (int): number of nodes in MIS
    """
    options = {
        "WLSACCESSID" : "55c88bbb-cfc0-4d09-b763-e37e0d3bd302",
        "WLSSECRET" : "f21e8941-0ff2-49c2-8500-6883437c1bb2",
        "LICENSEID" : 2677447,

    }
    env = Env(params=options)
    model = Model("arc_mis", env=env)

    x = model.addVars(node_dict, vtype=GRB.BINARY)

    model.addConstrs(x[i] + x[j] <= 1 for i, j in arcs)
    model.setObjective(quicksum(x[i] for i in node_dict), GRB.MAXIMIZE)
    model.setParam("OutputFlag", 0)
    model.setParam("TimeLimit", time_limit)
    model.optimize()

    result_vals = model.getAttr("x", x)
    optimal_nodes = [node_dict[i] for i in node_dict if round(result_vals[i]) == 1]
    mis_size = len(optimal_nodes)

    return optimal_nodes, mis_size, model.MIPGap, model.Runtime


def save_mis_csv(optimal_nodes, output_path):
    with open(output_path, 'w', newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(optimal_nodes)


def plot_3d(points):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for p in points:
        ax.scatter(p[0], p[1], p[2], s=10, c='blue')
    plt.show()


def plot_plotly(points):
    x, y, z = zip(*points)
    fig = px.scatter_3d(x=x, y=y, z=z)
    fig.show()


if __name__ == '__main__':
    # Default file inputs for testing
    node_file = "../data/node_data006.pkl"
    arc_file = "../data/arc_data006.pkl"
    output_csv = "../data/Optimal_Points_006.csv"

    with open(node_file, 'rb') as f:
        node_dict = pickle.load(f)

    with open(arc_file, 'rb') as f:
        arcs = pickle.load(f)

    optimal_nodes, mis_size = optimize_arc_mis(node_dict, arcs, time_limit=600, show_plot=True)

    save_mis_csv(optimal_nodes, output_csv)

    print(f"✔ MIS size: {mis_size}")
    print(f"✔ Saved {len(optimal_nodes)} optimal points to {output_csv}")

    # Optional plot
    plot_plotly(optimal_nodes)
