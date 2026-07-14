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
def dual(x_vals, nodes, arcs, Voxels, A, D, Beamlets, Struct_list, S_dict,
         max_constrs={"dummy": 1000000}, min_constrs={"dummy": -1000000}):
    model = Model('dual_problem')
    x_vals = np.array(x_vals)

    # Define variables using Gurobi variables
    x = model.addVars(nodes, vtype=GRB.CONTINUOUS, name="x")
    y = model.addVars(Beamlets, vtype=GRB.CONTINUOUS, name="y")
    z = model.addVars(Voxels, vtype=GRB.CONTINUOUS, name="z")

    # Set objective function to maximize the sum of z
    model.setObjective(-10*sum(x_vals) + 1000 + sum(z[i] for i in range(Voxels)), GRB.MINIMIZE)

    count = 0
    for key, val in nodes.items():
        model.addConstr(x[key] == x_vals[count])
        count += 1

    # Store delta constraints in a dict for easier tracking
    delta_constrs = {}
    for (i, j) in arcs:
        delta_constrs[(i, j)] = model.addConstr(x[i] + x[j] <= 1)

    # Constraints for pi_plus and pi_minus
    pi_plus_constrs = []
    pi_minus_constrs = []
    for v in range(Voxels):
        # Define the constraints for pi_plus and pi_minus
        Ay = sum(A[v,b] * y[b] for b in range(Beamlets))
        Dx = sum(D[v,n] * x[n] for n in nodes)
        pi_plus_constr = model.addConstr(z[v] >= Ay - Dx)    ##sum(A[v, k] * y[k] for k in range(Beamlets)) - sum(D[v, k] * x[k] for k in range(len(x_vals))))
        pi_minus_constr = model.addConstr(z[v] >= -Ay + Dx)  #-(sum(A[v, k] * y[k] for k in range(Beamlets)) - sum(D[v, k] * x[k] for k in range(len(x_vals)))))

        pi_plus_constrs.append(pi_plus_constr)
        pi_minus_constrs.append(pi_minus_constr)



    # Optimize the model
    model.setParam("OutputFlag", 0)
    model.optimize()

    if model.status == GRB.OPTIMAL:
        # Extract dual values in the format you requested
        pi_plus_duals = {f"pi_plus_{v}": pi_plus_constrs[v].Pi for v in range(Voxels)}
        pi_minus_duals = {f"pi_minus_{v}": pi_minus_constrs[v].Pi for v in range(Voxels)}
        delta_duals = {f"delta_{i}_{j}": [(i, j), delta_constrs[(i, j)].Pi] for (i, j) in arcs}

        with open('variable_values.csv', mode='w', newline='') as file:
            writer = csv.writer(file)

            for i, var in z.items():
                writer.writerow([i, var.X])

        # Return the dual variables in the format you requested
        return (
            pi_plus_duals,  # pi_plus duals
            pi_minus_duals,  # pi_minus duals
            delta_duals,  # delta duals
            model.ObjVal,  # Objective value
            "optimal"  # Status
        )
    elif model.status == GRB.UNBOUNDED:
        raise Exception("crap unbounded")
    elif model.status == GRB.INFEASIBLE:
        raise Exception("crap infeasible")
    else:
        raise Exception("Dual problem not solved correctly.")



def master(cuts, feas_cuts, nodes, arcs, voxels, D, pi_plus, pi_minus, delta,
           max_constrs, min_constrs, optimality_cut):
    model = Model("master problem")

    x = model.addVars(nodes, vtype=GRB.BINARY, name="x")
    theta = model.addVar(vtype=GRB.CONTINUOUS, name="theta")
    l = 10
    model.setObjective(-l * x.sum() + theta + 1000, GRB.MINIMIZE)

    # Pairwise constraints
    model.addConstrs(x[i] + x[j] <= 1 for (i, j) in arcs)
    val = 0
    val = LinExpr()

    count = 0

    if optimality_cut:
        # π terms (pi_plus and pi_minus)
        for v in range(voxels):
            val += (-pi_plus[f"pi_plus_{v}"] + pi_minus[f"pi_minus_{v}"]) * sum(D[v, i] * x[i] for i in nodes)
            count +=1


        # δ (delta constraints for adjacency)
        for s, stuff in delta.items():
            (i,j), dual_val = stuff
            val += dual_val * (x[i] + x[j] - 1)
            count += 1

        model.addConstr(theta >= val)

    else:
        raise Exception("No optimality cut made")

    # Add lazy optimality cuts (from previous iterations)
    for cut_data, const in cuts:
        cut_expr = LinExpr()
        for idx, coeff in cut_data:
            cut_expr += coeff * x[idx]
        model.addConstr(theta >= cut_expr + const)

    print(f"cut len: {len(cuts)}")

    if len(cuts) == 10:
        print("hi")

    model.addConstr(theta >= -1000000000)  # Avoid unbounded theta

    model.setParam("OutputFlag", 0)
    model.optimize()

    if model.status == GRB.OPTIMAL or model.status == GRB.TIME_LIMIT:
        mis = np.array([x[i].X for i in nodes])

        # Extract current optimality cut
        cut_data = []
        cut_constant = 0
        for i in range(val.size()):
            var = val.getVar(i)
            coeff = val.getCoeff(i)
            if var.VarName.startswith("x"):
                index = int(var.VarName.split("[")[1].split("]")[0])
                cut_data.append((index, coeff))
            else:
                cut_constant += coeff * var.X  # Evaluated constant

        return mis, model.objVal, (cut_data, cut_constant)

    if model.status == GRB.UNBOUNDED:
        raise Exception("Master problem is unbounded")



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
    for i in range(7):
        for j in range(7):
            nodes[(i+ 3)+ 10* (j+3)] = (0,0,0)
    arcs = []
    arc_matrix = create_arc_matrix()
    for i in range(len(arc_matrix)):
        for j in range(len(arc_matrix[j])):
            if arc_matrix[i][j] == 1 and i != j:
                arcs.extend([(i,j)])
    new_arcs = []
    for (i,j) in arcs:
        if i in nodes.keys() and j in nodes.keys():
            new_arcs.append((i,j))
    A = np.array(df)
    A = A.transpose()
    D = np.array(D)*10

    return A, D, nodes, new_arcs

A, D, nodes, arcs = load_toy_model()

roi_list = ["light", "plum", "tumor"]
voxel_dict = {"light" : [1,2,3,11,12,13,21,22,23,31,32,33,41,42,43,51,52,53,61,62,63], "plum" : [91,81,71,92,93,82,83,72,73],"tumor" : [88,87,79,78,77,76,75,69,68,67,66,65,64, 59,58,57,56,55, 49,48,47,46, 88, 89,94,95,96,97,98,99]}
max_constrs = {"light":7, "plum": 5}
min_constrs = {"light": 4, "tumor": 4}
avg_constrs = {"avg":100, "max":100}
dvh_constrs = {"dvh":[.3,3]}




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

pi_plus, pi_minus, delta, dual_obj, status = dual(x_vals=mis,arcs = arcs,nodes=nodes, Voxels=A.shape[0], A=A, D=D, Beamlets=A.shape[1], Struct_list=roi_list, S_dict=voxel_dict, max_constrs=max_constrs, min_constrs=min_constrs)
optimality_cut = True

if status == "unbounded":
    # Dual is unbounded, generate feasibility cut
    cut_data, cut_constant = build_feasibility_cut(pi_plus, pi_minus, mu_u, mu_l, mis, D, max_constrs, min_constrs)
    feasibility_cuts.append((cut_data, cut_constant))
    optimality_cut = False
    print(f"adding extreme ray cut")
elif status == "infeasible":
    print("Dual problem is infeasible, skipping feasibility cut generation.")

l = 10

while gap != TOL and iteration <= MAX_ITERS:
    iteration += 1

    # Call the master problem with the current cuts
    x_vals, master_obj, cut = master(
        cuts=cuts, feas_cuts=feasibility_cuts, nodes=nodes, arcs=arcs, voxels=A.shape[0], D=D,
        pi_plus=pi_plus, pi_minus=pi_minus, delta=delta,
        max_constrs=max_constrs, min_constrs=min_constrs, optimality_cut=optimality_cut
    )

    cuts.append(cut)

    # Call the dual problem to get the dual variables
    pi_plus, pi_minus, delta, dual_obj, status = dual(
        x_vals=x_vals, nodes = nodes, arcs=arcs, Voxels=A.shape[0], A=A, D=D, Beamlets=A.shape[1],
        Struct_list=roi_list, S_dict=voxel_dict, max_constrs=max_constrs, min_constrs=min_constrs
    )

    if status == "unbounded":
        # Dual is unbounded, generate feasibility cut
        cut_data, cut_constant = build_feasibility_cut(pi_plus, pi_minus,  x_vals, D, max_constrs, min_constrs)
        feasibility_cuts.append((cut_data, cut_constant))
        optimality_cut = False
        print(f"adding extreme ray cut")

    elif status == "infeasible":
        print("Dual problem is infeasible, skipping feasibility cut generation.")
        break  # Optionally, handle infeasibility gracefully (e.g., stop the loop)

    else:
        gap = -1* (master_obj - dual_obj) / (abs(master_obj) + 1e-6)
        optimality_cut = True
    print(f"Master Obj: {master_obj}")
    print(f"Dual obj: {dual_obj}")
    print(f"Gap: {gap:.6f}")
    print(f"")

print(f"Converged after {iteration} iterations with gap {gap:.6f}")
