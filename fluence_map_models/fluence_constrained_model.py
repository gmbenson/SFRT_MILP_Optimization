from gurobipy import *
import numpy as np
import pandas as pd
import pickle
import csv

def presolve(nodes, arcs, time_limit):
    options = {
        "WLSACCESSID": "4540b5e9-11c6-443c-aa3a-2e5b2dfcbb86",
        "WLSSECRET": "0203e490-6c19-4ed9-9cdc-b733f200e726",
        "LICENSEID": 2544834,

    }
    env = Env(params=options)
    model = Model("presolve", env=env)

    x = model.addVars(nodes, vtype=GRB.BINARY)

    model.setObjective(x.sum(), GRB.MAXIMIZE)
    model.addConstrs(x[i] + x[j] <= 1 for (i,j) in arcs)

    model.Params.TimeLimit = time_limit
    model.update()
    model.optimize()

    mis_size = 0
    for i in nodes.keys():  # Iterate over the
        if x[i].X >= 0.5:
            mis_size += 1
    print("Rough MIS set size calculated for a maximum of " + str(time_limit) + " seconds")
    print("MIS found with size " + str(mis_size))

    return mis_size

def get_bounds(bound, num_runs = 3, bound_search_ratio = 4):
    bounds = [bound - (bound/(i+bound_search_ratio)) for i in range(num_runs-1)]
    bounds.append(bound)
    return bounds


def initial_optimization(A, D, nodes, arcs, beamlets = 8682, voxels = 212367, mis_size = 3):
    #beamlets = 10
    options = {
        "WLSACCESSID" : "4540b5e9-11c6-443c-aa3a-2e5b2dfcbb86",
        "WLSSECRET" : "0203e490-6c19-4ed9-9cdc-b733f200e726",
        "LICENSEID" : 2544834,

    }
    env = Env(params=options)
    model = Model('fluence_map_model', env=env)

    y = model.addMVar(beamlets, vtype=GRB.CONTINUOUS)
    x = model.addMVar(len(nodes), vtype=GRB.BINARY)
    lam = 10
    peak_weight = 10
    D = peak_weight * D


    z = model.addVars(voxels, vtype=GRB.CONTINUOUS)

    for i in range(voxels):
        model.addConstr(z[i] >= (A[i,:] @ y - peak_weight* D[i,:] @ x))  # t[i] >= A@y - D@x
        #model.addConstr(z >= -(A @ y - D @ x))
        model.addConstr(z[i] >= -(A[i, :] @ y - peak_weight* D[i, :] @ x))  # t[i] >= -(A@y - D@x)

    # Set the objective function using quicksum
    model.setObjective(z.sum(), GRB.MINIMIZE)
    #model.setObjective(-1 * lam * x.sum() + y.transpose() @ A.transpose() @ A @ y - 2* D.transpose() @ x.transpose() @ A @ y + D.transpose()@x.transpose() @ D @ x, GRB.MINIMIZE)

    model.addConstr(x.sum() >= mis_size)

    model.addConstrs(x[i] + x[j] <= 1 for i,j in arcs)

    rhs = np.full(voxels, 30)

    for i in range(voxels):
        model.addConstr(A[i, :] @ y <= rhs[i])

    model.Params.TimeLimit = 43200
    model.update()
    model.optimize()

    dose = 0
    beamlet_use = 0
    deviation = 0

    if model.status == GRB.OPTIMAL:
        # Open a CSV file to write the values

        for i in range(beamlets):  # Iterate over the
            dose += y[i].X
            if i != 0:
                beamlet_use += 1
        #
        for i in range(voxels):  # Iterate over the
            deviation += z[i].X

    if dose == 0 and deviation == 0 and beamlet_use == 0:
        dose = "model infeasible"
        deviation = "model infeasible"

    return dose, deviation, beamlet_use

def main_optimization(A, D, nodes, arcs, beamlets = 8682, voxels = 212367, fluence_bound = 0):
    #beamlets = 10
    options = {
        "WLSACCESSID" : "4540b5e9-11c6-443c-aa3a-2e5b2dfcbb86",
        "WLSSECRET" : "0203e490-6c19-4ed9-9cdc-b733f200e726",
        "LICENSEID" : 2544834,

    }
    env = Env(params=options)
    model = Model('fluence_map_model', env=env)

    y = model.addMVar(beamlets, vtype=GRB.CONTINUOUS)
    x = model.addMVar(len(nodes), vtype=GRB.BINARY)
    lam = 10
    peak_weight = 10
    D = peak_weight * D

    z = model.addVars(voxels, vtype=GRB.CONTINUOUS)

    for i in range(voxels):
        model.addConstr(z[i] >= (A[i,:] @ y - peak_weight* D[i,:] @ x))  # t[i] >= A@y - D@x
        #model.addConstr(z >= -(A @ y - D @ x))
        model.addConstr(z[i] >= -(A[i, :] @ y - peak_weight* D[i, :] @ x))  # t[i] >= -(A@y - D@x)

    # Set the objective function using quicksum
    model.setObjective(x.sum(), GRB.MAXIMIZE)
    #model.setObjective(-1 * lam * x.sum() + y.transpose() @ A.transpose() @ A @ y - 2* D.transpose() @ x.transpose() @ A @ y + D.transpose()@x.transpose() @ D @ x, GRB.MINIMIZE)

    model.addConstr(z.sum() <= fluence_bound)

    model.addConstrs(x[i] + x[j] <= 1 for i,j in arcs)

    rhs = np.full(voxels, 30)

    for i in range(voxels):
        model.addConstr(A[i, :] @ y <= rhs[i])

    model.Params.TimeLimit = 43200
    model.update()
    model.optimize()

    dose = 0
    beamlet_use = 0
    deviation = 0
    node_count = 0
    if model.status == GRB.OPTIMAL:
        # Open a CSV file to write the values
        for i in range(len(nodes)):  # Iterate over the
            if x[i].X >= 0.5:
                node_count += 1
        for i in range(beamlets):  # Iterate over the
            dose += y[i].X
            #if y[i] != 0:
            #    beamlet_use += 1
        #
        for i in range(voxels):  # Iterate over the
            deviation += z[i].X

    if dose == 0 and deviation == 0 and beamlet_use == 0:
        dose = "model infeasible"
        deviation = "model infeasible"

    with open("fluence_constr_results.csv", 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(dose)
        writer.writerow(deviation)
        writer.writerow(node_count)

    return node_count, dose, deviation, beamlet_use
def get_adjacent_squares(x, y, grid_size=10):
    adjacent_squares = []

    # Directions: up, down, left, right, and diagonals
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1),  # vertical and horizontal
                   (0,0)]  # diagonal

    # Check all possible directions
    for dx, dy in directions:
        new_x, new_y = x + dx, y + dy

        # Check if the new position is within grid bounds
        if 0 <= new_x < grid_size and 0 <= new_y < grid_size:
            adjacent_squares.append((new_x, new_y))

    return adjacent_squares

def get_adjacent_nodes(x, y, grid_size=10):
    adjacent_squares = []

    # Directions: up, down, left, right, and diagonals
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1),  # vertical and horizontal
                  (1,1),(-1,1),(1,-1),(-1,-1),(1,2),
                  (-1,2), (1,-2), (-1,-2), (2,1), (-2,1),
                  (2,-1), (-2,-1), (2,2),(2,-2), (-2,2), (-2,-2),
                  (2,0), (-2,0), (0,2), (0,-2)]  # diagonal

    # Check all possible directions
    for dx, dy in directions:
        new_x, new_y = x + dx, y + dy

        # Check if the new position is within grid bounds
        if 0 <= new_x < grid_size and 0 <= new_y < grid_size:
            adjacent_squares.append((new_x, new_y))

    return adjacent_squares
def create_arc_matrix(grid_size=10):
    # Initialize a grid with zeros (100x100 matrix)
    D = np.zeros((grid_size * grid_size, grid_size * grid_size), dtype=int)

    # Iterate over each square in the grid
    for x in range(grid_size):
        for y in range(grid_size):
            # Calculate the index for the current square (row and column in 1D array)
            current_index = x * grid_size + y
            # Get the adjacent squares
            adjacent = get_adjacent_nodes(x, y, grid_size)

            # Mark the current square and adjacent squares with 1
            D[current_index, current_index] = 1  # A square is adjacent to itself
            for ax, ay in adjacent:
                adjacent_index = ax * grid_size + ay
                D[current_index, adjacent_index] = 1  # Mark as adjacent

    return D
def create_adjacency_matrix(grid_size=10):
    # Initialize a grid with zeros (100x100 matrix)
    D = np.zeros((grid_size * grid_size, grid_size * grid_size), dtype=int)

    # Iterate over each square in the grid
    for x in range(grid_size):
        for y in range(grid_size):
            # Calculate the index for the current square (row and column in 1D array)
            current_index = x * grid_size + y
            # Get the adjacent squares
            adjacent = get_adjacent_squares(x, y, grid_size)

            # Mark the current square and adjacent squares with 1
            D[current_index, current_index] = 1  # A square is adjacent to itself
            for ax, ay in adjacent:
                adjacent_index = ax * grid_size + ay
                D[current_index, adjacent_index] = 1  # Mark as adjacent

    return D

def load_toy_model():
    df = pd.read_excel(r"C:\Users\Grant\Downloads\Toy_Model.xlsx")
    df = df.transpose()
    df = df.values.tolist()
    voxels = range(len(df[0]))
    A = []
    for i in range(len(df)):
        holder = []
        for j in range(len(df[i])):
            holder.append(df[i][j])

    D = create_adjacency_matrix()

    nodes = {}
    for i in range(len(df[0])):
        nodes[i] = (0,0,0)
    arcs = []
    arc_matrix = create_arc_matrix()
    for i in range(len(arc_matrix)):
        for j in range(len(arc_matrix[j])):
            if arc_matrix[i][j] == 1 and i != j:
                arcs.extend([(i,j)])
    A = np.array(df)
    A = A.transpose()
    D = np.array(D)

    return A, D, nodes, arcs


with open("../fluence_map_models/data/sparse_matrix.pkl", 'rb') as openfile:
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
A = A[:-1]
#A = inf_matrix.toarray()

#A, D, nodes, arcs = load_toy_model()

mis_size = presolve(nodes,arcs, 1200)
dose, deviaiton, beamlet_use = initial_optimization(A,D,nodes, arcs, voxels=A.shape[0],beamlets = A.shape[1], mis_size=mis_size)
results = []
bounds = get_bounds(deviaiton, 7, bound_search_ratio=2)
for bound in bounds:
    results.append(main_optimization(A, D, nodes, arcs, voxels = A.shape[0], beamlets = A.shape[1], fluence_bound = bound))#212001)


print("Optimization finished")
for result in results:
    print("mis size: " + str(result[0]) + ", total dose: " + str(result[1]) + ", total deviaiton: " + str(result[2]))