import pickle
import numpy as np
from gurobipy import *
import csv
import pandas as pd
from scipy.sparse import csr_matrix

def optimization(A, D, nodes, arcs, beamlets = 8682, voxels = 212367):
    beamlets = 10
    '''
    options = {
        "WLSACCESSID": "15bc6745-8f1e-4f24-9910-e232914fe388",
        "WLSSECRET": "1d2cd7b3-9c64-42cb-a67f-6410e887d6ca",
        "LICENSEID": 2544834,
    }
    env = Env(params=options)
    '''
    model = Model('fluence_map_model') #, env=env)

    y = model.addMVar(beamlets, vtype=GRB.CONTINUOUS)
    x = model.addMVar(len(nodes), vtype=GRB.BINARY)
    lam = 10
    peak_weight = 1
    D = peak_weight*D


    z = model.addVars(voxels, vtype=GRB.CONTINUOUS)

    for i in range(voxels):
        model.addConstr(z[i] >= (A[i,:] @ y - peak_weight* D[i,:] @ x))  # t[i] >= A@y - D@x
        #model.addConstr(z >= -(A @ y - D @ x))
        model.addConstr(z[i] >= -(A[i, :] @ y - peak_weight* D[i, :] @ x))  # t[i] >= -(A@y - D@x)

    # Set the objective function using quicksum
    model.setObjective(-1* lam * x.sum() + z.sum(), GRB.MINIMIZE)
    #model.setObjective(-1 * lam * x.sum() + y.transpose() @ A.transpose() @ A @ y - 2* D.transpose() @ x.transpose() @ A @ y + D.transpose()@x.transpose() @ D @ x, GRB.MINIMIZE)

    model.addConstrs(x[i] + x[j] <= 1 for i,j in arcs)

    rhs = np.full(voxels, 30)

    for i in range(voxels):
        model.addConstr(A[i, :] @ y <= rhs[i])

    model.Params.TimeLimit = 43200
    model.update()
    model.optimize()

    holder = []
    holder = model.getAttr('x')
    used_nodes = [holder[i] for i in range(len(nodes)) if
                  holder[i] > 0.5]  # Use a threshold of 0.5 to account for numerical precision

    print("used nodes: " + str(used_nodes))
    print("total num used: " + str(len(used_nodes)))

    if model.status == GRB.OPTIMAL:
        # Open a CSV file to write the values

        print("x vals") # Iterate over the
        print([x[i].X for i in range(voxels)])

        print("y vals")
        for i in range(beamlets):  # Iterate over the
            print(([y[i].X]))
        #
        print("z vals")
        for i in range(voxels):  # Iterate over the
            print(([z[i].X]))

        with open('nodes_solution.csv', mode='w', newline='') as file:
            writer = csv.writer(file)
            # Loop through the MVars (x, y, z) and write the values to the CSV
            for i in range(len(nodes)):  # Iterate over the
                writer.writerow([x[i].X])

        with open('fluence_solution.csv', mode='w', newline='') as file:
            writer = csv.writer(file)
            # Loop through the MVars (x, y, z) and write the values to the CSV
            for i in range(beamlets):  # Iterate over the rows
                writer.writerow([y[i].X])
        with open('z_solution.csv', mode='w', newline='') as file:
            writer = csv.writer(file)
            # Loop through the MVars (x, y, z) and write the values to the CSV
            exit()
            for i in range(voxels):  # Iterate over the rows
                writer.writerow([z[i].X])


        print("Solution saved to 'solution.csv'")
    else:
        print("Optimization was not successful.")


    return holder, used_nodes


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
'''
with open("../fluence_map_models/data/sparse_matrix.pkl", 'rb') as openfile:
    A = pickle.load(openfile)
print("Fluence matrix loaded")
with open("../fluence_map_models/data/sparse_target_matrix.pkl", 'rb') as openfile:
    D = pickle.load(openfile)
'''
'''
nonzero = D.getnnz(axis=1)
for i, nonzeroes in enumerate(nonzero):
    if nonzeroes != 0:
        print(f"Row {i} has {nonzeroes} non-zero entries")
exit()
'''
'''
print("Target matrix loaded")
with open("../fluence_map_models/data/nodes.pkl", 'rb') as openfile:
    nodes = pickle.load(openfile)
with open("../fluence_map_models/data/arcs.pkl", 'rb') as openfile:
    arcs = pickle.load(openfile)
print("Nodes and arcs loaded")
A = A[:212001]
#D = D[123500:132001]
#A = inf_matrix.toarray()
'''
A, D, nodes, arcs = load_toy_model()

solution, nodes = optimization(A, D, nodes, arcs, voxels = 100)#212001)

print("Optimization finished")

with open('fluence_sols.csv', 'w', newline="") as csvfile:
    writer = csv.writer(csvfile, delimiter=',')
    writer.writerow(solution)

print(len(solution))