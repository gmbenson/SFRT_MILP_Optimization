from gurobipy import *
import numpy as np
import pandas as pd
import pickle


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

def dual(x, Voxels, A, D, Beamlets, Struct_list, S_dict, max_constrs = {"dummy" : 1000000}, min_constrs = {"dummy": -1000000}):

    # Create the optimization model
    model = Model('dual_problem')
    Voxels = [i for i in range(Voxels)]
    Beamlets = [i for i in range(Beamlets)]
    x = np.array(x)

    #Struct_list.append("dummy")
    #S_dict["dummy"] = [0]


    # Create variables
    pi_plus = model.addVars(Voxels, vtype=GRB.CONTINUOUS, name="pi_plus", lb=0)
    pi_minus = model.addVars(Voxels, vtype=GRB.CONTINUOUS, name="pi_minus", lb=0)
    mu_u = model.addVars(max_constrs, vtype=GRB.CONTINUOUS, name="mu_u", lb=0)
    mu_l = model.addVars(min_constrs, vtype=GRB.CONTINUOUS, name="mu_l", lb=0)

    # Objective function
    obj = LinExpr()

    # Summation over v for pi_plus and pi_minus
    for v in Voxels:
        obj += -pi_plus[v] * (D[v, :] @ x) + pi_minus[v] * (D[v, :] @ x)

    # Summation over s for mu, nu, gamma, delta
    for u, maxi in max_constrs.items():
        obj += -mu_u[u] * maxi

    for l, mini in min_constrs.items():
        obj += mu_l[l] * mini

    model.setObjective(obj, GRB.MAXIMIZE)

    # Constraints
    # Constraint (1): pi_plus + pi_minus <= 1 for all v in V
    for v in Voxels:
        model.addConstr(pi_plus[v] + pi_minus[v] <= 1)

    # Constraint (2): The sum over v and s
    for b in Beamlets:
        lhs = LinExpr()
        for v in Voxels:
            lhs += (-A[v, b] * pi_plus[v] + A[v, b] * pi_minus[v])
        for s in max_constrs:
            for v in S_dict[s]:
                lhs += A[v, b] * -mu_u[s]
        for s in min_constrs:
            for v in S_dict[s]:
                lhs += A[v, b] * mu_l[s]
        model.addConstr(lhs <= 0)

    model.setParam("Outputflag", 0)
    model.optimize()
    if model.status == GRB.OPTIMAL:
        return (
            {k: pi_plus[k].X for k in pi_plus},
            {k: pi_minus[k].X for k in pi_minus},
            {k: mu_u[k].X for k in mu_u},
            {k: mu_l[k].X for k in mu_l},
            model.ObjVal,
            "optimal"
        )
    elif model.status == GRB.UNBOUNDED:
        print("Dual is unbounded — extracting extreme ray for feasibility cut...")

        # Extract unbounded ray
        ray_pi_plus = {k: pi_plus[k].UnbdRay for k in pi_plus}
        ray_pi_minus = {k: pi_minus[k].UnbdRay for k in pi_minus}
        ray_mu_u = {k: mu_u[k].UnbdRay for k in mu_u}
        ray_mu_l = {k: mu_l[k].UnbdRay for k in mu_l}

        return ray_pi_plus, ray_pi_minus, ray_mu_u, ray_mu_l, None, "unbounded"

    else:
        raise Exception("dual problem not solved correctly")
def build_feasibility_cut(ray_pi_plus, ray_pi_minus, ray_mu_u, ray_mu_l, x_vals, D, max_constrs, min_constrs):
    cut_expr = LinExpr()  # Initialize the linear expression for the cut
    rhs_violation = 0  # This will be the constant part of the cut

    # Loop over all x variables (nodes)
    for i in range(len(x_vals)):
        coeff = 0
        # Loop over all voxels (rows in D matrix)
        for v in range(D.shape[0]):
            coeff += (-ray_pi_plus[v] + ray_pi_minus[v]) * D[v, i]  # Interaction of dual variables with x[i] via D

        # Contributions from max constraints
        for s in max_constrs:
            if s in ray_mu_u:  # Ensure the constraint exists in the dictionary
                for v in voxel_dict[s]:  # Look up max constraints based on voxel dictionary
                    coeff += -ray_mu_u[s] * D[v, i]  # Max constraint contribution

        # Contributions from min constraints
        for s in min_constrs:
            if s in ray_mu_l:  # Ensure the constraint exists in the dictionary
                for v in voxel_dict[s]:  # Look up min constraints based on voxel dictionary
                    coeff += ray_mu_l[s] * D[v, i]  # Min constraint contribution

        cut_expr += coeff * x_vals[i]  # Add to the cut expression

    # Right-hand side violation (constant part of the cut)
    rhs_violation = 0
    for v in range(D.shape[0]):
        dose = sum(D[v, i] * x_vals[i] for i in range(len(x_vals)))  # Compute dose for voxel v
        rhs_violation += (-ray_pi_plus[v] + ray_pi_minus[v]) * dose  # Contributions from dual variables

    for s, max_val in max_constrs.items():
        if s in ray_mu_u:
            rhs_violation += -ray_mu_u[s] * max_val  # Contributions from max constraints

    for s, min_val in min_constrs.items():
        if s in ray_mu_l:
            rhs_violation += ray_mu_l[s] * min_val  # Contributions from min constraints

    return cut_expr, rhs_violation  # Return the cut expression and its constant



def master(cuts, feas_cuts, nodes, arcs, voxels, D, pi_plus, pi_minus, mu_u, mu_l, max_constrs, min_constrs, optimality_cut): #= {"dummy" : 1000000}, min_constrs = {"dummy": -1000000}):
    model = Model("master problem")

    voxels = [i for i in range(voxels)]

    x = model.addVars(len(nodes), vtype=GRB.BINARY, name = "x")
    theta = model.addVar( vtype=GRB.CONTINUOUS, name = "theta")
    l = 10
    model.setObjective(- l * x.sum() + theta, GRB.MINIMIZE)

    model.addConstrs(x[i] + x[j] <= 1 for (i,j) in arcs)
    if optimality_cut:
        val = LinExpr()

        # Summation over v for pi_plus and pi_minus
        for v in voxels:
            for i in range(len(nodes)):
                val += (-pi_plus[v] + pi_minus[v]) * D[v, i] * x[i]

        # Summation over s for mu, nu, gamma, delta
        for u, maxi in max_constrs.items():
            val += -mu_u[u] * maxi

        for l, mini in min_constrs.items():
            val += mu_l[l] * mini
        model.addConstr(theta >= val)
    else:
        val = LinExpr()
        raise Exception("no optimality cut made")

    if cuts:
        for cut_data, const in cuts:
            cut_expr = LinExpr()
            for idx, coeff in cut_data:
                cut_expr += coeff * x[idx]
            model.addConstr(theta >= cut_expr + const)

    if feas_cuts:
        for cut_expr, const in feas_cuts:
            model.addConstr(0 >= cut_expr + const)

    print(f"cut len: {len(cuts)}")

    model.addConstr(theta >= -1000000000)


    model.setParam("OutputFlag", 0)
    model.optimize()

    if model.status == GRB.OPTIMAL or model.status == GRB.TIME_LIMIT:
        mis = [x[i].X for i in range(len(nodes))]
        mis = np.array(mis)

        cut_data = []
        cut_constant = 0

        # Safe access: iterate over all terms in val
        for i in range(val.size()):
            var = val.getVar(i)
            coeff = val.getCoeff(i)

            # Only keep x-variable terms in the cut
            if var.VarName.startswith("x"):
                index = int(var.VarName.split('[')[1].split(']')[0])  # x[5] → 5
                cut_data.append((index, coeff))
            else:
                cut_constant += coeff * var.X  # fixed value from dual variables

        return mis, model.objVal, (cut_data, cut_constant)
    if model.status == GRB.UNBOUNDED:
        raise Exception("master unbounded")



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
    df = df.values.tolist()[:10]
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

A, D, nodes, arcs = load_toy_model()
roi_list = ["max", "avg", "dvh"]
voxel_dict = {"max" : [1,2,3,4,5,11,12,13,14], "avg" : [9,10,19,20,70],"dvh" : [70,71,72,73,74,80,81,82,83,98,97,96,95,99]}
max_constrs = {"max":100}
min_constrs = {"dvh": 5}
avg_constrs = {"avg":100, "max":100}
dvh_constrs = {"dvh":[.3,3]}



TOL = 1e-4  # convergence tolerance
MAX_ITERS = 50
iteration = 0
feasibility_cuts = []  # List to store feasibility cuts]
optimality_cuts = []
TOL = 1e-4  # convergence tolerance
MAX_ITERS = 50
iteration = 0

cuts = []
# Main loop to solve the master and dual problems
gap = float('inf')
old_obj = float('-inf')
mis_size, mis = presolve(nodes,arcs, 60)

feasibility_cuts = []

pi_plus, pi_minus, mu_u, mu_l, dual_obj, status = dual(x=mis, Voxels=A.shape[0], A=A, D=D, Beamlets=A.shape[1], Struct_list=roi_list, S_dict=voxel_dict, max_constrs=max_constrs, min_constrs=min_constrs)
optimality_cut = True

if status == "unbounded":
    # Dual is unbounded, generate feasibility cut
    cut_data, cut_constant = build_feasibility_cut(pi_plus, pi_minus, mu_u, mu_l, mis, D, max_constrs, min_constrs)
    feasibility_cuts.append((cut_data, cut_constant))
    optimality_cut = False
    print(f"adding extreme ray cut")
elif status == "infeasible":
    print("Dual problem is infeasible, skipping feasibility cut generation.")



while gap != TOL and iteration <= MAX_ITERS:
    iteration += 1

    # Call the master problem with the current cuts
    x_vals, master_obj, cut = master(
        cuts=cuts, feas_cuts= feasibility_cuts,
        nodes=nodes, arcs=arcs, voxels=A.shape[0], D=D,
        pi_plus=pi_plus, pi_minus=pi_minus, mu_u=mu_u, mu_l=mu_l,
        max_constrs=max_constrs, min_constrs=min_constrs, optimality_cut = optimality_cut
    )

    cuts.append(cut)

    # Call the dual problem to get the dual variables
    pi_plus, pi_minus, mu_u, mu_l, dual_obj, status = dual(
        x=x_vals, Voxels=A.shape[0], A=A, D=D, Beamlets=A.shape[1],
        Struct_list=roi_list, S_dict=voxel_dict,
        max_constrs=max_constrs, min_constrs=min_constrs
    )

    if status == "unbounded":
        # Dual is unbounded, generate feasibility cut
        cut_data, cut_constant = build_feasibility_cut(pi_plus, pi_minus, mu_u, mu_l, x_vals, D, max_constrs, min_constrs)
        feasibility_cuts.append((cut_data, cut_constant))
        optimality_cut = False
        print(f"adding extreme ray cut")

    elif status == "infeasible":
        print("Dual problem is infeasible, skipping feasibility cut generation.")
        break  # Optionally, handle infeasibility gracefully (e.g., stop the loop)

    else:
        # Feasible dual solution, compute gap
        gap = (master_obj - dual_obj) / (abs(master_obj) + 1e-6)
        optimality_cut = True
    print(f"Master Obj: {master_obj}")
    print(f"Dual obj: {dual_obj}")
    print(f"Gap: {gap:.6f}")
    print(f"")

print(f"Converged after {iteration} iterations with gap {gap:.6f}")
