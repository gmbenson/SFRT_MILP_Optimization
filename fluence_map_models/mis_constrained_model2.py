from gurobipy import *
import numpy as np
import pandas as pd
import pickle
import csv

def presolve(nodes, arcs, time_limit):
    options = {
        "WLSACCESSID" : "4540b5e9-11c6-443c-aa3a-2e5b2dfcbb86",
        "WLSSECRET" : "0203e490-6c19-4ed9-9cdc-b733f200e726",
        "LICENSEID" : 2544834,
    }
    #env = Env(params=options)
    model = Model("presolve")#, env=env)
    x = model.addVars(nodes, vtype=GRB.BINARY, name = "x")

    model.setObjective(x.sum(), GRB.MAXIMIZE)
    model.addConstrs(x[i] + x[j] <= 1 for (i,j) in arcs)

    model.Params.TimeLimit = time_limit
    model.update()
    model.optimize()
    if model.status == GRB.OPTIMAL or model.status == GRB.TIME_LIMIT:
        mis_size = 0
        mis = []
        for i in x.keys():  # Iterate over the
            mis.append(x[i].X)
            if x[i].X >= 0.5:
                mis_size += 1
        print("Rough MIS set size calculated for a maximum of " + str(time_limit) + " seconds")
        print("MIS found with size " + str(mis_size))

    return mis_size, mis

def find_mis_sets(mis_size, num_runs = 3):
    mis_sets = [mis_size - i for i in range(num_runs)]
    return mis_sets


def optimization(A, D, nodes, arcs, beamlets = 8682, voxels = 212367, mis_size = 3, mis = [], roi_list = [], voxel_list = [], max_constrs = {}, avg_constrs = {}, dvh_constrs = {}):
    #beamlets = 10
    options = {
        "WLSACCESSID" : "4540b5e9-11c6-443c-aa3a-2e5b2dfcbb86",
        "WLSSECRET" : "0203e490-6c19-4ed9-9cdc-b733f200e726",
        "LICENSEID" : 2544834,

    }
    #env = Env(params=options)
    model = Model('fluence_map_model')#, env=env)

    y = model.addMVar(beamlets, vtype=GRB.CONTINUOUS)
    x = model.addMVar(len(nodes), vtype=GRB.BINARY)
    lam = 1
    peak_weight = (66/15)
    D = peak_weight * D
    z = model.addMVar(voxels, vtype=GRB.CONTINUOUS)

    print("vars loaded")

    model.addConstr(z >= A@y - peak_weight*D@x)
    model.addConstr(z >= -1*(A @ y - peak_weight * D @ x))
    '''
    z = model.addVars(voxels, vtype=GRB.CONTINUOUS)    
   
    for i in range(voxels):
        model.addConstr(z[i] >= (A[i,:] @ y - peak_weight* D[i,:] @ x))  # t[i] >= A@y - D@x
        #model.addConstr(z >= -(A @ y - D @ x))
        model.addConstr(z[i] >= -(A[i, :] @ y - peak_weight* D[i, :] @ x))  # t[i] >= -(A@y - D@x)
    '''
    # Set the objective function using quicksum
    #model.setObjective(z.sum(), GRB.MINIMIZE)
    model.setObjective(z.sum(), GRB.MINIMIZE)

    #model.setObjective(-1 * lam * x.sum() + y.transpose() @ A.transpose() @ A @ y - 2* D.transpose() @ x.transpose() @ A @ y + D.transpose()@x.transpose() @ D @ x, GRB.MINIMIZE)

    model.addConstr(x.sum() >= mis_size)

    model.addConstrs(x[i] + x[j] <= 1 for i,j in arcs)

    print("basic model loaded")

    for struct, max_val in max_constrs.items():
        for voxel in voxel_list[roi_list.index(struct)]:
            model.addConstr((A@y)[voxel] <= max_val)
    print("max constriaints loaded")

    for struct, av in avg_constrs.items():
        expr = LinExpr()
        for voxel in voxel_list[roi_list.index(struct)]:
            expr += (A@y)[voxel]
        model.addConstr(expr/len(voxel_list[roi_list.index(struct)]) <= av)

    print("average constraints loaded")

    zeta = {}
    dummy = {}
    for struct, vals in dvh_constrs.items():
        ratio, dose = vals
        u_vox = voxel_list[roi_list.index(struct)][-1]
        l_vox = voxel_list[roi_list.index(struct)][0]
        zeta[struct] = model.addVar(name=struct, ub= np.infty, lb=-1*np.infty, vtype=GRB.CONTINUOUS)
        dummy[struct] = model.addMVar(u_vox-l_vox, vtype=GRB.CONTINUOUS)
        model.addConstr(dummy[struct] >= (A@y)[l_vox:u_vox] - zeta[struct])
        model.addConstr(dummy[struct] >= 0)

        model.addConstr(zeta[struct]+(1/((1-ratio)*(u_vox-l_vox)))*dummy[struct].sum() <= dose)

    print("dvh constraints loaded")

    model.Params.TimeLimit = 7200
    model.update()
    mis = np.array(mis)


    #for i in range(len(mis)):
    x.Start = mis

    model.optimize()

    dose = 0
    beamlet_use = 0
    deviation = 0
    print("Y")
    print(y.X)


    #if model.status == GRB.OPTIMAL:
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
        beamlet_use = "model infeasible"
    with open("mis_constr_dose_results.csv", 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([dose])
    with open("mis_constr_dev_results.csv", 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([deviation])
    with open("mis_constr_beamlet_results.csv", 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([beamlet_use])

    print(f"total deviation: {sum(deviation)}")

    return dose, deviation, beamlet_use


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
    D = np.array(D)*10

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
#A = A[:212000]

data = []
with open("../fluence_map_models/data/inf_params.pkl", 'rb') as openfile:
    data = pickle.load(openfile)

voxel_pos, voxel_dat, tumor_voxel, all_voxels = data
'''
A, D, nodes, arcs = load_toy_model()
roi_list = ["max", "avg", "dvh"]
voxel_list = [[1,2,3,4,5,11,12,13,14],  [9,10,19,20,70], [70,71,72,73,74,80,81,82,83,98,97,96,95,99]]
max_constrs = {"max":100}
avg_constrs = {"avg":100, "max":100}
dvh_constrs = {"dvh":[.3,3]}
'''
roi_list = ["GTV", "LUNGS_NOT_GTV", "PTV", "ESOPHAGUS", "HEART", "LUNG_L", "LUNG_R", "CORD", "SKIN", "BODY"]
voxel_list = all_voxels
max_constrs = {"GTV":69, "PTV":69, "ESOPHAGUS":66, "HEART":27, "LUNG_L":66, "LUNG_R":66, "CORD":50, "SKIN":60,"LUNGS_NOT_GTV":66}
avg_constrs = {"ESOPHAGUS": 34, "LUNGS_NOT_GTV": 21}
dvh_constrs = {"ESOPHAGUS": [.17,60], "HEART":[.5,30], "LUNGS_NOT_GTV": [.37,20]}



#mis_sets = find_mis_sets(presolve(nodes,arcs, 60), num_runs=8)
mis_size, mis = presolve(nodes,arcs, 60)
results = []
optimization(A,D,nodes,arcs, voxels=A.shape[0], beamlets=A.shape[1], mis_size=mis_size, mis = mis, roi_list=roi_list, voxel_list=voxel_list, max_constrs=max_constrs, avg_constrs=avg_constrs, dvh_constrs=dvh_constrs)
#for mis_set in mis_sets:
#    results.append([mis_set, optimization(A, D, nodes, arcs, voxels = A.shape[0], beamlets = A.shape[1], mis_size= mis_set)])#212001)


print("Optimization finished")
for result in results:
    print("mis min: " + str(result[0]) + ", total dose: " + str(result[1][0]) + ", total deviaiton: " + str(result[1][1]) + ", num beamlets used: " + str(result[1][2]))