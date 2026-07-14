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

def dual(x, Voxels, A, D, Beamlets, Struct_list, S_dict, max_constrs = {"dummy" : 1000000}, min_constrs = {"dummy": -1000000}, u_mean_constrs = {"dummy":100000000}, l_mean_constrs = {"dummy":-1000000}, u_dvh_constrs = {"dummy":[.5, 10000000]}, l_dvh_constrs = {"dummy": [.5,-100000]}):

    # Create the optimization model
    model = Model('dual_problem')
    Voxels = [i for i in range(Voxels)]
    Beamlets = [i for i in range(Beamlets)]
    x = np.array(x)

    Struct_list.append("dummy")
    S_dict["dummy"] = [0]

    u_omega_dict = [(s,v) for s in u_dvh_constrs for v in S_dict[s]]
    l_omega_dict = [(s, v) for s in l_dvh_constrs for v in S_dict[s]]

    # Create variables
    pi_plus = model.addVars(Voxels, vtype=GRB.CONTINUOUS, name="pi_plus", lb=0)
    pi_minus = model.addVars(Voxels, vtype=GRB.CONTINUOUS, name="pi_minus", lb=0)
    mu_u = model.addVars(max_constrs, vtype=GRB.CONTINUOUS, name="mu_u", lb=0)
    mu_l = model.addVars(min_constrs, vtype=GRB.CONTINUOUS, name="mu_l", lb=0)
    delta_u = model.addVars(u_mean_constrs, vtype=GRB.CONTINUOUS, name="delta_u", lb=0)
    delta_l = model.addVars(l_mean_constrs, vtype=GRB.CONTINUOUS, name="delta_l", lb=0)
    sigma_u = model.addVars(u_dvh_constrs, vtype=GRB.CONTINUOUS, name="sigma_U", lb=0)
    sigma_l = model.addVars(l_dvh_constrs, vtype=GRB.CONTINUOUS, name="sigma_L", lb=0)
    omega_tau = model.addVars(u_omega_dict, vtype=GRB.CONTINUOUS, name="omega_tau", lb=0)
    omega_bar_tau = model.addVars(l_omega_dict, vtype=GRB.CONTINUOUS, name="omega_bar_tau", lb=0)

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

    for u, maxi in u_mean_constrs.items():
        obj += -delta_u[u] * maxi

    for l, mini in l_mean_constrs.items():
        obj += delta_l[l] * mini

    for u, maxi in u_dvh_constrs.items():
        obj += -sigma_u[u] * maxi[1]

    for l, mini in l_dvh_constrs.items():
        obj += sigma_l[l] * mini[1]

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
        for s in u_mean_constrs:
            for v in S_dict[s]:
                lhs += A[v, b] * (1 / len(S_dict[s])) * delta_u[s]
        for s in l_mean_constrs:
            for v in S_dict[s]:
                lhs += A[v, b] * -(1 / len(S_dict[s])) * delta_l[s]
        for s,v in u_omega_dict:
            lhs += A[v, b] * -omega_tau[s,v]
        for s,v in l_omega_dict:
            for v in S_dict[s]:
                lhs += A[v, b] * omega_bar_tau[s,v]
        model.addConstr(lhs <= 0)


    # Constraint (3) and (4): Summations over alpha in A^s for sigma_U and sigma_L
    for s, [alpha, val] in u_dvh_constrs.items():
        for v in S_dict[s]:
            model.addConstr(omega_tau[s, v] - (1 / ((1 - alpha) * len(S_dict[s]))) * sigma_u[s] <= 0)

    for s, [alpha, val] in l_dvh_constrs.items():
        for v in S_dict[s]:
            model.addConstr(omega_bar_tau[s, v] - (1 / ((1 - alpha) * len(S_dict[s]))) * sigma_l[s] <= 0)

    # Constraint (5): Relation between sigma_U and omega_bar_tau
    for s in u_dvh_constrs:
        model.addConstr(-sigma_u[s] + sum(omega_tau[s, v] for v in S_dict[s]) == 0)

    # Constraint (6): Relation between sigma_L and omega_tau
    for s  in l_dvh_constrs:
        model.addConstr(sigma_l[s] + sum(omega_bar_tau[s, v] for v in S_dict[s]) == 0)

    # Optimize the model

    model.setParam("Outputflag", 0)
    model.optimize()

    return (
        {k: pi_plus[k].X for k in pi_plus},
        {k: pi_minus[k].X for k in pi_minus},
        {k: mu_u[k].X for k in mu_u},
        {k: mu_l[k].X for k in mu_l},
        {k: delta_u[k].X for k in delta_u},
        {k: delta_l[k].X for k in delta_l},
        {k: sigma_u[k].X for k in sigma_u},
        {k: sigma_l[k].X for k in sigma_l},
        {k: omega_tau[k].X for k in omega_tau},
        {k: omega_bar_tau[k].X for k in omega_bar_tau},
        model.ObjVal
    )



def master(cuts, nodes, arcs, voxels, D, pi_plus, pi_minus, mu_u, mu_l, delta_u, delta_l, sigma_u, sigma_l, max_constrs = {"dummy" : 1000000}, min_constrs = {"dummy": -1000000}, u_mean_constrs = {"dummy":100000000}, l_mean_constrs = {"dummy":-1000000}, u_dvh_constrs = {"dummy":[.5, 10000000]}, l_dvh_constrs = {"dummy": [.5,-100000]}):
    model = Model("master problem")

    voxels = [i for i in range(voxels)]

    x = model.addVars(len(nodes), vtype=GRB.BINARY, name = "x")
    theta = model.addVar( vtype=GRB.CONTINUOUS, name = "theta")
    l = 10
    model.setObjective(l * x.sum() - theta, GRB.MAXIMIZE)

    model.addConstrs(x[i] + x[j] <= 1 for (i,j) in arcs)

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

    for u, maxi in u_mean_constrs.items():
        val += -delta_u[u] * maxi

    for l, mini in l_mean_constrs.items():
        val += delta_l[l] * mini

    for u, maxi in u_dvh_constrs.items():
        val += -sigma_u[u] * maxi[1]

    for l, mini in l_dvh_constrs.items():
        val += sigma_l[l] * mini[1]

    if cuts:
        for cut_data, const in cuts:
            cut_expr = LinExpr()
            for idx, coeff in cut_data:
                cut_expr += coeff * x[idx]
            model.addConstr(theta >= cut_expr + const)

    model.addConstr(theta >= val)



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
    else:
        raise RuntimeError("Master problem did not solve to optimality.")



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

A, D, nodes, arcs = load_toy_model()
roi_list = ["max", "avg", "dvh"]
voxel_dict = {"max" : [1,2,3,4,5,11,12,13,14], "avg" : [9,10,19,20,70],"dvh" : [70,71,72,73,74,80,81,82,83,98,97,96,95,99]}
max_constrs = {"max":100}
avg_constrs = {"avg":100, "max":100}
dvh_constrs = {"dvh":[.3,3]}

mis_size, mis = presolve(nodes,arcs, 60)

pi_plus, pi_minus, mu_u, mu_l, delta_u, delta_l, sigma_u, sigma_l, omega_tau, omega_bar_tau, s_obj = dual(x=mis, Voxels=A.shape[0], A=A, D=D, Beamlets=A.shape[1], Struct_list=roi_list, S_dict=voxel_dict, max_constrs=max_constrs, u_mean_constrs=avg_constrs, u_dvh_constrs=dvh_constrs)

TOL = 1e-4  # convergence tolerance
MAX_ITERS = 50
iteration = 0

cuts = []
gap = float('inf')
old_obj = float('-inf')

while gap != TOL and iteration <= 1000:
    iteration+=1

    x_vals, master_obj, cut = master(cuts = cuts,
        nodes=nodes, arcs=arcs, voxels=A.shape[0], D=D,
        pi_plus=pi_plus, pi_minus=pi_minus, mu_u=mu_u, mu_l=mu_l,
        delta_u=delta_u, delta_l=delta_l, sigma_u=sigma_u, sigma_l=sigma_l,
        max_constrs=max_constrs, u_mean_constrs=avg_constrs, u_dvh_constrs=dvh_constrs
    )

    cuts.append(cut)


    pi_plus, pi_minus, mu_u, mu_l, delta_u, delta_l, sigma_u, sigma_l, omega_tau, omega_bar_tau, dual_obj = dual(
        x=x_vals, Voxels=A.shape[0], A=A, D=D, Beamlets=A.shape[1],
        Struct_list=roi_list, S_dict=voxel_dict,
        max_constrs=max_constrs, u_mean_constrs=avg_constrs, u_dvh_constrs=dvh_constrs
    )

    gap = (master_obj - dual_obj) / (abs(master_obj) + 1e-6)

    print(f"Master Objective: {master_obj:.6f}")
    print(f"Dual Objective:   {dual_obj:.6f}")
    print(f"")
    print(f"Duality Gap:      {gap:.6f}")
    print(f"")


print(f"\nConverged after {iteration} iterations with gap {gap:.6f}")
print("Final MIS size:", int(sum(x_vals)))

print("final mis size: " + str(sum(x_vals)))