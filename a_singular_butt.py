from gurobipy import *
import numpy as np
import pandas as pd
import pickle
import csv


def presolve(nodes, arcs, time_limit):
    options = {
        "WLSACCESSID" : "55c88bbb-cfc0-4d09-b763-e37e0d3bd302",
        "WLSSECRET" : "f21e8941-0ff2-49c2-8500-6883437c1bb2",
        "LICENSEID" : 2677447,

    }
    env = Env(params=options)
    model = Model("presolve", env=env)
    x = model.addVars(nodes, vtype=GRB.BINARY, name = "x")

    model.setObjective(x.sum(), GRB.MAXIMIZE)
    model.addConstrs(x[i] + x[j] <= 1 for (i,j) in arcs)

    model.Params.TimeLimit = time_limit
    model.update()
    model.setParam("OutputFlag", 0)
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
def opt(x_vals, nodes, arcs, Voxels, A, D, Beamlets, Struct_list, S_dict,
         max_constrs={"dummy": 1000000}, min_constrs={"dummy": -1000000}, avg_constrs={"dummy": 1000000}, dvh_constrs ={"dummy": -1000000}):
    options = {
        "WLSACCESSID": "55c88bbb-cfc0-4d09-b763-e37e0d3bd302",
        "WLSSECRET": "f21e8941-0ff2-49c2-8500-6883437c1bb2",
        "LICENSEID": 2677447,

    }
    env = Env(params=options)
    model = Model("opt", env=env)
    x_vals = np.array(x_vals)

    # Define variables using Gurobi variables
    x = model.addVars(nodes, vtype=GRB.BINARY, name="x")
    y = model.addVars(Beamlets, vtype=GRB.CONTINUOUS, name="y")
    z = model.addVars(Voxels, vtype=GRB.CONTINUOUS, name="z")

    # Set objective function to maximize the sum of z
    model.setObjective(-150*sum(x[i] for i in nodes) + sum(z[i] for i in range(Voxels)), GRB.MINIMIZE)

    model.addConstr(sum(x[i] for i in nodes) == 4)

    # Store delta constraints in a dict for easier tracking
    delta_constrs = {}
    for (i, j) in arcs:
        delta_constrs[(i, j)] = model.addConstr(x[i] + x[j] <= 1)

    # Constraints for pi_plus and pi_minus
    pi_plus_constrs = []
    pi_minus_constrs = []

    scaling = 10

    struct_targets = {"funny": 2, "light": 5, "plum": 1, "tumor": 6, "border":6}

    for v in range(Voxels):
        target = 0
        for st, val in struct_targets.items():
            if v in S_dict[st]:
                target = val
        # Define the constraints for pi_plus and pi_minus
        Ay = sum(A[v,b] * y[b] for b in range(Beamlets))
        Dx = sum(D[v,n] * x[n] for n in nodes)
        pi_plus_constr = model.addConstr(z[v] >= Ay - 3*Dx - target)##sum(A[v, k] * y[k] for k in range(Beamlets)) - sum(D[v, k] * x[k] for k in range(len(x_vals))))
        pi_minus_constr = model.addConstr(z[v] >= -Ay + 3*Dx + target) #-(sum(A[v, k] * y[k] for k in range(Beamlets)) - sum(D[v, k] * x[k] for k in range(len(x_vals)))))

        pi_plus_constrs.append(pi_plus_constr)
        pi_minus_constrs.append(pi_minus_constr)


    mu_u_constrs = {}
    for struct, max_val in max_constrs.items():
        for i in S_dict[struct]:
            # Define constraints for mu_u
            expr = sum(A[i, k] * y[k] for k in range(Beamlets))
            c = model.addConstr(expr <= max_val)

    # Constraints for mu_l (lower bounds on A @ y)
    mu_l_constrs = {}
    for struct, min_val in min_constrs.items():
        for i in S_dict[struct]:
            # Define constraints for mu_l
            expr = sum(A[i, k] * y[k] for k in range(Beamlets))
            c = model.addConstr(expr >= min_val)

    for struct, avg_val in avg_constrs.items():
        expr = LinExpr()
        for i in S_dict[struct]:
            expr += sum(A[i, k] * y[k] for k in range(Beamlets))
        model.addConstr(expr/len(S_dict[struct]) <= avg_val)

    zeta = {}
    dummy = {}

    for struct, (ratio, dose) in dvh_constrs.items():
        # Continuous slack variable
        zeta[struct] = model.addVar(name=f"zeta_{struct}", lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)

        # Dummy variables for each voxel in structure
        dummy_vars = []
        for i in S_dict[struct]:
            # A@y at voxel i
            dose_expr = sum(A[i, k] * y[k] for k in range(Beamlets))

            # Dummy variable represents max(0, dose_expr - zeta)
            dummy_i = model.addVar(lb=0, vtype=GRB.CONTINUOUS, name=f"dummy_{struct}_{i}")
            model.addConstr(dummy_i >= dose_expr - zeta[struct])
            dummy_vars.append(dummy_i)

        # DVH constraint
        model.addConstr(
            zeta[struct] + (1 / ((1 - ratio) * (len(S_dict[struct])))) * quicksum(dummy_vars) <= dose,
            name=f"dvh_constr_{struct}"
        )

    # Optimize the model
    model.setParam("OutputFlag", 0)
    model.optimize()

    if model.status == GRB.OPTIMAL:
        # Extract dual values in the format you requested
        #pi_plus_duals = {f"pi_plus_{v}": pi_plus_constrs[v].Pi for v in range(Voxels)}
        #pi_minus_duals = {f"pi_minus_{v}": pi_minus_constrs[v].Pi for v in range(Voxels)}
        #delta_duals = {f"delta_{i}_{j}": [(i, j), delta_constrs[(i, j)].Pi] for (i, j) in arcs}
        with open('deviation_values.csv', mode='w', newline='') as file:
            writer = csv.writer(file)

            for i, var in z.items():
                writer.writerow([i, var.X])

        with open('dose_values.csv', mode='w', newline='') as file:
            writer = csv.writer(file)
            v_num = 0
            for i, var in z.items():
                writer.writerow([v_num, sum(A[i, k] * y[k].X for k in range(Beamlets))])
                v_num +=1

        print(f"num nodes:")
        print(sum(x[i].X for i in nodes))
        print(f"")
        print(f"total deviation:")
        print(sum(z[i].X for i in range(Voxels)))
        print(f"beamlet intensity:")
        print([y[i].X for i in range(Beamlets)])
        print(f"total dose:")
        print(sum(A[i, k] * y[k].X for k in range(Beamlets) for i in range(Voxels)))

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
def create_arcs_from_nodes(node_set, grid_size=10):
    arcs = set()
    node_set = set(node_set)

    # 8-connectivity: up/down/left/right + diagonals
    directions = [
        (dx, dy)
        for dx in range(-2, 3)
        for dy in range(-2, 3)
        if not (dx == 0 and dy == 0)
    ]

    for node in node_set:
        x, y = node % grid_size, node // grid_size

        for dx, dy in directions:
            nx, ny = x + dx, y + dy
            if 0 <= nx < grid_size and 0 <= ny < grid_size:
                neighbor_index = ny * grid_size + nx
                if neighbor_index in node_set:
                    arcs.add(tuple(sorted((node, neighbor_index))))  # avoid (i,j) and (j,i)

    return list(arcs)


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

    D = create_adjacency_matrix()
    arc_matrix = create_arc_matrix()

    # Step 1: Define original tumor indices (before shrinking)
    original_tumor = {}
    for i in range(7):
        for j in range(6):
            idx = (i + 3) + 10 * (j + 3)
            original_tumor[idx] = (0, 0, 0)
    for remove_idx in [33, 34, 43]:
        original_tumor.pop(remove_idx, None)

    tumor_set = set(original_tumor)
    border = set()

    # Step 2: Shrink tumor
    for voxel in tumor_set:
        x, y = voxel % 10, voxel // 10
        neighbors = get_adjacent_squares(x, y)
        for nx, ny in neighbors:
            neighbor_index = ny * 10 + nx
            if neighbor_index not in tumor_set:
                border.add(voxel)
                break

    shrunken_tumor_set = tumor_set - border
    tumor_border = list(border)

    # Step 3: Create nodes and arcs only for the shrunken tumor
    nodes = {i: (0, 0, 0) for i in [37, 38, 45, 46, 47, 48, 54, 55, 56, 57, 58, 64, 65, 66, 67, 68, 74, 75, 76, 77, 78, 85,86]}
    arcs = create_arcs_from_nodes(nodes)

    # Format dose and adjacency matrix
    A = np.array(df).T
    D = np.array(D) * 10

    return A, D, nodes, arcs, list(shrunken_tumor_set), tumor_border


A, D, nodes, arcs, tumor, tumor_border = load_toy_model()

roi_list = ["light", "plum", "funny", "tumor", "border"]



roi_list = ["light", "plum", "funny", "tumor", "border"]
voxel_dict = {"funny":[90,80,70,60,50,40,30], "light" : [1,2,3,11,12,13,21,22,23,31,32,33,41,42,43,51,52,53,61,62,63], "plum" : [91,81,71,92,93,82,83,72,73],'tumor': [37, 38, 45, 46, 47, 48, 54, 55, 56, 57, 58, 64, 65, 66, 67, 68, 74, 75, 76, 77, 78, 85,86], 'border': [27,28, 49, 59,69,79, 35, 36, 39, 73, 44, 83, 84, 53, 95, 96, 87, 89, 88, 63]}
print(voxel_dict)
max_constrs = {"light":10, "plum": 11, "funny":12}
min_constrs = {"light": 0, "tumor": 0}
avg_constrs = {"light":4, "plum":5}
dvh_constrs = {"light":[.3,6], "plum":[.3,4]}




feasibility_cuts = []  # List to store feasibility cuts]
optimality_cuts = []
TOL = 1e-4  # convergence tolerance
MAX_ITERS = 300
iteration = 0

cuts = []
# Main loop to solve the master and dual problems
gap = float('inf')
old_obj = float('-inf')
mis_size, mis = presolve(nodes,arcs, 60)

feasibility_cuts = []

opt(x_vals=mis,arcs = arcs,nodes=nodes, Voxels=A.shape[0], A=A, D=D, Beamlets=A.shape[1], Struct_list=roi_list, S_dict=voxel_dict, max_constrs=max_constrs, min_constrs=min_constrs, avg_constrs=avg_constrs, dvh_constrs=dvh_constrs)

exit()
