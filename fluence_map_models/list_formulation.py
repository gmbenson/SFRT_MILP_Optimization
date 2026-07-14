import pickle
import numpy as np
from gurobipy import *
import csv

def optimization(A, D, nodes, arcs):
    options = {
        "WLSACCESSID": "15bc6745-8f1e-4f24-9910-e232914fe388",
        "WLSSECRET": "1d2cd7b3-9c64-42cb-a67f-6410e887d6ca",
        "LICENSEID": 2544834,
    }
    #env = Env(params=options)
    model = Model('fluence_map_model')#, env=env)

    D = np.array(D)

    max_val = 0
    for inner_dict in A.values():
        if inner_dict:
            max_val = max(max_val, max(inner_dict.keys()))

    y = model.addVars(range(max_val+1), vtype=GRB.CONTINUOUS)
    x = model.addVars(nodes, vtype=GRB.BINARY)
    z = model.addVars(len(A), vtype= GRB.BINARY)


    print("variables created")

    lam = 100000
    target_weight = 1000

    for v in range(len(A)):
        for b in A[v].keys():
            for n in nodes.keys():
                if n in D and v in D[n]:
                    model.addConstr(z[v] >= A[v][b] * y[b] - target_weight*D[n][v] * x[n])
                else:
                    model.addConstr(z[v] >= A[v][b] * y[b])

    model.setObjective(-lam * quicksum(x[n] for n in nodes) + quicksum(z[v]**2 for v in range(len(A))), GRB.MINIMIZE)

    #model.setObjective(lam * quicksum(x[n] for n in nodes) + quicksum(
    #    (A[v][b] * y[b] - target_weight*D[n][v] * x[n]) ** 2 for v in range(len(A)) for b in A[v].keys()  for n in range(len(nodes))), GRB.MINIMIZE)

    print("Objective function set")

    model.addConstrs(x[i] + x[j] <= 1 for i,j in arcs)
    model.addConstrs(quicksum(A[v][b] * y[b] for b in range(A[v].keys())) <= 30 for v in range(len(A)))

    print("constraints set")

    model.Params.TimeLimit = 43200
    model.update()
    model.optimize()

    holder = []
    holder = model.getAttr('x')
    return holder



with open("../fluence_map_models/data/inf_list.pkl", 'rb') as openfile:
    A = pickle.load(openfile)
print("Fluence matrix loaded")
with open("../fluence_map_models/data/target_list.pkl", 'rb') as openfile:
    D = pickle.load(openfile)
print("Target matrix loaded")
with open("../fluence_map_models/data/nodes.pkl", 'rb') as openfile:
    nodes = pickle.load(openfile)
with open("../fluence_map_models/data/arcs.pkl", 'rb') as openfile:
    arcs = pickle.load(openfile)
print("Nodes and arcs loaded")

#A = inf_matrix.toarray()

solution = optimization(A, D, nodes, arcs)

print("Optimization finished")

with open('fluence_sols.csv', 'w', newline="") as csvfile:
    writer = csv.writer(csvfile, delimiter=',')
    writer.writerows(solution)

print(len(solution))