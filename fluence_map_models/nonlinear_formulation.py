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

    y = model.addVars(len(A[0]), vtype=GRB.CONTINUOUS)
    x = model.addVars(nodes, vtype=GRB.BINARY)

    lam = 1000

    model.setObjective(lam * quicksum(x[n] for n in nodes) + quicksum(
        (A[v, b] * y[b] - D[n,v] * x[n]) ** 2 for b in range(len(A[0])) for v in range(len(A)) for n in range(len(nodes))), GRB.MINIMIZE)

    model.addConstrs(x[i] + x[j] <= 1 for i,j in arcs)
    model.addConstrs(quicksum(A[v, b] * y[b] for b in range(len(A[0]))) <= 30 for v in range(len(A)))

    model.Params.TimeLimit = 43200
    model.update()
    model.optimize()

    holder = []
    holder = model.getAttr('x')
    return holder



with open("../fluence_map_models/data/inf_matrix.pkl", 'rb') as openfile:
    A = pickle.load(openfile)
print("Fluence matrix loaded")
with open("../fluence_map_models/data/sparse_target_matrix.pkl", 'rb') as openfile:
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