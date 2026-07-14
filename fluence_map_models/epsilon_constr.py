import pickle
import numpy as np
from gurobipy import *

def node_optimization(A, D, nodes, arcs, fluence_bound, time_limit):
    model = Model('fluence_map_model')

    A = np.array(A)
    D = np.array(D)

    y = model.addVars(len(A[0]), vtype=GRB.CONTINUOUS)
    x = model.addVars(nodes, vtype=GRB.BINARY)

    model.setObjective(quicksum(x[n] for n in nodes), GRB.MAXIMIZE)

    model.addConstrs(x[i] + x[j] <= 1 for i,j in arcs)
    model.addConstrs(quicksum(A[i, j] * y[j] for j in range(len(A[0]))) <= 30 for i in range(len(A)))
    model.addConstrs(quicksum((A[v, b] * y[b] - D[v,n] * x[n]) ** 2 for b in range(len(A[0])) for v in range(len(A)) for n in range(len(nodes))) <= fluence_bound)

    model.Params.TimeLimit = time_limit
    model.update()
    model.optimize()

    x_vals = {}
    y_vals = []

    for x in model.getVars():
        if x.VarName.startswith("x"):
            x_vals[x.VarName] = x.X  # Use getAttr to get the value
        elif x.VarName.startswith("y"):
            y_vals.append(x.X)
        else:
            raise Exception("NOOOOO")
    return x_vals, y_vals, model.ObjVal

def fluence_optimization(A, D, nodes, arcs, num_nodes, time_limit):
    model = Model('fluence_map_model')

    A = np.array(A)
    D = np.array(D)

    y = model.addVars(len(A[0]), vtype=GRB.CONTINUOUS)
    x = model.addVars(nodes, vtype=GRB.BINARY)

    lam = 1000

    model.setObjective(quicksum((A[v, b] * y[b] - D[v,n] * x[n]) ** 2 for b in range(len(A[0])) for v in range(len(A)) for n in range(len(nodes))), GRB.MINIMIZE)

    model.addConstrs(x[i] + x[j] <= 1 for i,j in arcs)
    model.addConstrs(quicksum(A[i, j] * y[j] for j in range(len(A[0]))) <= 30 for i in range(len(A)))
    model.addConstr(quicksum(x[n] for n in nodes) >= num_nodes)
    model.Params.TimeLimit = time_limit
    model.update()
    model.optimize()

    x_vals = {}
    y_vals = []

    for x in model.getVars():
        if x.VarName.startswith("x"):
            x_vals[x.VarName] = x.X  # Use getAttr to get the value
        elif x.VarName.startswith("y"):
            y_vals.append(x.X)
        else:
            raise Exception("NOOOOO")
    return x_vals, y_vals, model.ObjVal

def run_model(A,D,nodes, arcs, num_iterations, iter_time = 600):
    num_nodes = 0
    fluence_bound = 0
    for i in range(num_iterations):
        x_vals, y_vals, fluence_bound = fluence_optimization(A,D,nodes,arcs,num_nodes, time_limit=iter_time)
        x_vals, y_vals, num_nodes = node_optimization(A,D,nodes,arcs, fluence_bound, time_limit=iter_time)
    return x_vals, y_vals




with open("../fluence_map_models/data/inf_matrix.pkl", 'rb') as openfile:
    A = pickle.load(openfile)
with open("../fluence_map_models/data/target_matrix.pkl", 'rb') as openfile:
    D = pickle.load(openfile)
with open("../fluence_map_models/data/nodes.pkl", 'rb') as openfile:
    nodes = pickle.load(openfile)
with open("../fluence_map_models/data/arcs.pkl", 'rb') as openfile:
    arcs = pickle.load(openfile)

#A = inf_matrix.toarray()

opt_nodes, opt_fluences = run_model(A, D, nodes, arcs, iter_time=600)
print("found " + str(len(opt_nodes)) + " nodes")
print(opt_nodes)