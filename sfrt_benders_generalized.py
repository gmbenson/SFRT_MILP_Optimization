import gurobipy as gp
from gurobipy import GRB
from gurobipy import *
import numpy as np
import pandas as pd
import pickle
import csv
from typing import List, Tuple

import gurobipy as gp
from gurobipy import GRB
import numpy as np

def benders_binary(A, B, b, c, f, arcs, max_iters=5000, tol=1e-6):
    """
    Benders Decomposition for min c^T x + f^T y s.t. A x + B y >= b,
    with binary x (complicating) and continuous y (easy).

    Parameters:
        A, B, b, c, f: problem data (numpy arrays)
        y_bounds: list of (lb, ub) for each y variable
        max_iters: max number of iterations
        tol: convergence tolerance

    Returns:
        Optimal x and y, and objective value
    """
    m, nx = A.shape
    _, ny = B.shape

    # Master Problem: Binary x, continuous theta
    master = gp.Model("Benders_Master_BinaryX")
    x = master.addVars(nx, vtype=GRB.BINARY, name="x")
    theta = master.addVar(lb=-GRB.INFINITY, name="theta")
    master.setObjective(gp.quicksum(c[i]*x[i] for i in range(nx)) + theta, GRB.MINIMIZE)
    master.addConstr(theta >= -10000000)
    master.addConstrs(x[i]+x[j] <=1 for (i,j) in arcs)
    master.setParam("OutputFlag", 0)

    for it in range(max_iters):
        master.optimize()
        if master.Status != GRB.OPTIMAL:
            print("Master not optimal.")
            break

        x_val = np.array([x[i].X for i in range(nx)])

        # Subproblem with fixed x
        sub = gp.Model("Benders_Subproblem")
        y = sub.addVars(ny, lb=0, name="y", vtype=GRB.CONTINUOUS)

        rhs = b - A @ x_val
        constraints = [sub.addConstr(gp.quicksum(B[i, j] * y[j] for j in range(ny)) >= rhs[i])
                       for i in range(m)]
        sub.setObjective(gp.quicksum(f[j] * y[j] for j in range(ny)), GRB.MINIMIZE)
        sub.setParam("OutputFlag", 0)

        # Solve the subproblem (dual variables will be stored in .Pi)
        sub.optimize()

        if sub.Status == GRB.OPTIMAL:
            y_val = np.array([y[j].X for j in range(ny)])
            sub_obj = f @ y_val
            theta_val = theta.X

            # Dual variables corresponding to the subproblem constraints
            duals = np.array([constraint.Pi for constraint in constraints])
            # Lagrange multipliers for constraints
            #print(f"Duals: {duals}")
            master_obj_val = c @ x_val + theta_val
            gap = abs(sub_obj - theta_val) / sub_obj
            if gap > tol:
                print(f"Iter {it}: master obj value = {master_obj_val}, theta = {theta_val:.4f}, sub_obj = {sub_obj:.4f}, gap = {abs(sub_obj - theta_val)/sub_obj:.6f}")
                # Add Benders optimality cut
                cut_expr = sum(duals[i] * (b[i] - gp.quicksum(A[i, k] * x[k] for k in range(nx)))
                               for i in range(m))
                master.addConstr(theta >= cut_expr)
            else:
                print(f"Converged in {it+1} iterations.")
                return x_val, y_val, c @ x_val + f @ y_val

        elif sub.Status == GRB.INFEASIBLE:
            print("Subproblem infeasible, adding feasibility cut.")
            sub.computeIIS()
            sub.write("subproblem.ilp")  # Optional: to inspect infeasibility

            duals = np.array([-c.FarkasDual for c in constraints])
            cut_expr = sum(duals[i] * (b[i] - gp.quicksum(A[i, k] * x[k] for k in range(nx)))
                           for i in range(m))
            master.addConstr(cut_expr <= 0)
        else:
            print("Subproblem not solved properly.")
            break

    print("Benders did not converge.")
    return None, None, None


def extract_benders_components(
    model: Model,
    x_vars: List[Var],
    all_vars: List[Var]
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Extracts components (c, f, A, B, b) for Benders decomposition from a Gurobi model.
    Assumes:
        - x_vars are complicating binary variables
        - all_vars contains x_vars + all continuous (easy) variables like y, z, etc.

    Returns:
        c: Objective coefficients for x (shape: [nx])
        f: Objective coefficients for y/z/others (shape: [ny])
        A: Constraint coefficients for x (shape: [m, nx])
        B: Constraint coefficients for y/z/others (shape: [m, ny])
        b: Right-hand side vector (shape: [m])
    """
    model.update()

    # Partition easy variables (not in x_vars)
    x_varnames = set(v.VarName for v in x_vars)
    easy_vars = [v for v in all_vars if v.VarName not in x_varnames]

    nx = len(x_vars)
    ny = len(easy_vars)

    # Create index maps
    x_idx = {v.VarName: i for i, v in enumerate(x_vars)}
    y_idx = {v.VarName: i for i, v in enumerate(easy_vars)}

    # Objective coefficients
    c = np.array([v.Obj for v in x_vars])
    f = np.array([v.Obj for v in easy_vars])

    A_rows = []
    B_rows = []
    b_vector = []

    for constr in model.getConstrs():
        expr = model.getRow(constr)
        sense = constr.Sense
        rhs = constr.RHS

        a_row = np.zeros(nx)
        b_row = np.zeros(ny)

        for i in range(expr.size()):
            var = expr.getVar(i)
            coeff = expr.getCoeff(i)

            if var.VarName in x_idx:
                a_row[x_idx[var.VarName]] = coeff
            elif var.VarName in y_idx:
                b_row[y_idx[var.VarName]] = coeff

        # Convert all constraints to >= form
        if sense == '=':
            A_rows.append(a_row)
            B_rows.append(b_row)
            b_vector.append(rhs)

            A_rows.append(-a_row)
            B_rows.append(-b_row)
            b_vector.append(-rhs)

        elif sense == '<':
            A_rows.append(-a_row)
            B_rows.append(-b_row)
            b_vector.append(-rhs)

        elif sense == '>':
            A_rows.append(a_row)
            B_rows.append(b_row)
            b_vector.append(rhs)

    A = np.vstack(A_rows)
    B = np.vstack(B_rows)
    b = np.array(b_vector)

    return c, f, A, B, b




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


def get_model(x_vals, nodes, arcs, Voxels, A, D, Beamlets, Struct_list, S_dict,
         max_constrs={"dummy": 1000000}, min_constrs={"dummy": -1000000}):
    model = Model('dual_problem')

    x_vals = np.array(x_vals)

    # Define variables using Gurobi variables
    x = model.addVars(nodes, vtype=GRB.BINARY, name="x")
    y = model.addVars(Beamlets, lb=0, vtype=GRB.CONTINUOUS, name="y")
    z = model.addVars(Voxels, lb=0, vtype=GRB.CONTINUOUS, name="z")

    # Set objective function to maximize the sum of z
    model.setObjective(-100*sum(x[i] for i in nodes) + sum(z[i] for i in range(Voxels)), GRB.MINIMIZE)

    delta_constrs = {}
    for (i, j) in arcs:
        delta_constrs[(i, j)] = model.addConstr(x[i] + x[j] <= 1)


    pi_plus_constrs = []
    pi_minus_constrs = []

    scaling = 10

    struct_targets = {"funny": 2}#, "light": 5, "plum": 1, "tumor": 10}

    for v in range(Voxels):
        target = 0
        for st, val in struct_targets.items():
            if v in S_dict[st]:
                target = val
        # Define the constraints for pi_plus and pi_minus
        Ay = sum(A[v,b] * y[b] for b in range(Beamlets))
        Dx = sum(D[v,n] * x[n] for n in nodes)
        pi_plus_constr = model.addConstr(z[v] >= Ay - 2*Dx - target)##sum(A[v, k] * y[k] for k in range(Beamlets)) - sum(D[v, k] * x[k] for k in range(len(x_vals))))
        pi_minus_constr = model.addConstr(z[v] >= -Ay + 2*Dx + target) #-(sum(A[v, k] * y[k] for k in range(Beamlets)) - sum(D[v, k] * x[k] for k in range(len(x_vals)))))

        pi_plus_constrs.append(pi_plus_constr)
        pi_minus_constrs.append(pi_minus_constr)

    model.update()

    return model




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
    for i in range(6):
        for j in range(6):
            nodes[(i+ 3)+ 10* (j+3)] = (0,0,0)
    nodes.pop(33)
    nodes.pop(43)
    nodes.pop(34)
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

roi_list = [ "funny"]#, "tumor","light", "plum"]
voxel_dict = {"funny":[0,1,2,3,4,5]}#[90,80,70,60,50,40,30]}#, #"light" : [1,2,3,11,12,13,21,22,23,31,32,33,41,42,43,51,52,53,61,62,63], "plum" : [91,81,71,92,93,82,83,72,73],"tumor" : [i for i in nodes]}
print(voxel_dict)
max_constrs = {"light":10, "plum": 6, "funny":12}
max_constrs = {"funny": 10}
min_constrs = {"funny": 1}
#min_constrs = {"light": 2, "tumor": 4}
avg_constrs = {"avg":100, "max":100}
dvh_constrs = {"dvh":[.3,3]}




mis_size, mis = presolve(nodes,arcs, 60)

model = get_model(x_vals=mis,arcs = arcs,nodes=nodes, Voxels=A.shape[0], A=A, D=D, Beamlets=A.shape[1], Struct_list=roi_list, S_dict=voxel_dict, max_constrs=max_constrs, min_constrs=min_constrs)


#x_vars = [model.getVarByName(f"x[{i}]") for i in nodes]
# Create consistent index mapping between node keys and variable positions
node_keys = list(nodes.keys())  # list of actual node IDs (e.g., 33, 34, etc.)
x_vars = [model.getVarByName(f"x[{i}]") for i in node_keys]
node_to_index = {node: idx for idx, node in enumerate(node_keys)}

y_vars = [model.getVarByName(f"y[{j}]") for j in range(A.shape[1])]
z_vars = [model.getVarByName(f"z[{v}]") for v in range(A.shape[0])]

all_vars = x_vars + y_vars + z_vars

c, f, A, B, b = extract_benders_components(model, x_vars, all_vars)

# Map arcs to indices compatible with x_vars
indexed_arcs = [(node_to_index[i], node_to_index[j]) for (i, j) in arcs if i in node_to_index and j in node_to_index]

benders_binary(A, B, b, c, f, indexed_arcs)
