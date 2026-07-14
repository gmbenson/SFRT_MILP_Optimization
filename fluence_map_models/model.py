import pickle
from gurobipy import *

def optimization(voxels, voxels_by_roi,tumor_voxels, arcs, neighborhood, beamlets, deposition,a,b):

    options = {
        "WLSACCESSID": "9a6728f5-f620-4829-b3b7-5e01b8e667d9",
        "WLSSECRET": "6f73eab7-cea6-4f6f-a167-510fd4aa9aa9",
        "LICENSEID": 2544834,
    }
    #env = Env(params=options)
    model = Model("optimal")#,env=env)  # create model

    m = 1000

    x = model.addVars(beamlets, vtype=GRB.CONTINUOUS, name="x")
    y = model.addVars(tumor_voxels, vtype=GRB.BINARY, name="y")
    p = model.addVars(voxels, vtype=GRB.CONTINUOUS, name="p")
    z_t = model.addVars(voxels, vtype=GRB.CONTINUOUS,name="z_t")
    z_o = model.addVars(voxels, vtype=GRB.CONTINUOUS,name="z_o")
    #d[v] = sum(deposition[(b,v)]*x[b] for b in beamlets)


    model.addConstrs(y[i] + y[j] <= 1 for i,j in arcs)

    model.addConstrs(y[i] + sum(y[j] for j in neighborhood[i]) >= 1 for i in tumor_voxels)

    model.addConstrs(
        p[v]/1000 >= a[k]*z_t[v] + b[k]*y[n]
        for n in tumor_voxels
        for v in {key: voxels[key] for key in neighborhood[n] if key in neighborhood[n]}.keys()
        for k in range(len(a))
    )

    model.addConstrs(
        z_t[v] <= m * y[n]
        for n in tumor_voxels
        for v in {key: voxels[key] for key in neighborhood[n] if key in neighborhood[n]}.keys()
    )
    model.addConstrs(z_t[v] <= sum(deposition[(b,v)]*x[b] for b in beamlets) for n in tumor_voxels for v in voxels.keys())
    model.addConstrs(
        z_t[v] >= sum(deposition[(b,v)]*x[b] for b in beamlets) - m*(1-y[n])
        for n in tumor_voxels
        for v in {key: voxels[key] for key in neighborhood[n] if key in neighborhood[n]}.keys()
    )




    model.addConstrs(
        p[v] >= 1 * z_o[v]
        for n in tumor_voxels
        for v in {key: voxels[key] for key in neighborhood[n] if key in neighborhood[n]}.keys()
        for k in range(len(a))
    )
    model.addConstrs(
        z_o[v] <= m * (1-y[n])
        for n in tumor_voxels
        for v in {key: voxels[key] for key in neighborhood[n] if key in neighborhood[n]}.keys()
    )
    model.addConstrs(
        z_o[v] <= sum(deposition[(b, v)] * x[b] for b in beamlets) for n in tumor_voxels for v in voxels.keys())

    model.addConstrs(
        z_o[v] >= sum(deposition[(b, v)] * x[b] for b in beamlets) - m * y[n]
        for n in tumor_voxels
        for v in {key: voxels[key] for key in neighborhood[n] if key in neighborhood[n]}.keys()
    )




    model.setObjective(quicksum(p[v] for v in voxels.keys()), GRB.MINIMIZE)
    model.setParam("OutputFlag", 1)
    model.Params.TimeLimit = 21600
    #model.Params.TimeLimit = 160
    model.update()
    model.optimize()
    if model.status == GRB.INFEASIBLE:
        # Compute the IIS (Irreducible Infeasible Subsystem)
        model.computeIIS()
        # Output the IIS to a file
        model.write("infeasible_model.ilp")
        print("IIS written to infeasible_model.ilp")
    model.write("model.mps")

    # Optionally, increase the level of detail in the log
    model.setParam('OutputFlag', 1)  # Enable detailed output during optimization
    model.setParam('LogFile', 'gurobi_log.txt')
    x_values = {}
    y_values = {}
    p_values = {}
    z_t_values = {}
    z_o_values = {}
    # Retrieve optimized values and split them into respective dictionaries
    for v in model.getVars():
        if v.VarName.startswith("x["):
            x_values[v.VarName] = v.X  # Use getAttr to get the value
        elif v.VarName.startswith("y["):
            y_values[v.VarName] = v.X  # Use getAttr to get the value
        elif v.VarName.startswith("p["):
            p_values[v.VarName] = v.X  # Store value in p_values
        elif v.VarName.startswith("z_t["):
            z_t_values[v.VarName] = v.X  # Use getAttr to get the value
        elif v.VarName.startswith("z_o["):
            z_o_values[v.VarName] = v.X  # Store value in p_values

    # Return the dictionaries for each variable type
    return x_values, y_values, p_values,z_t_values,z_o_values


def optimization2(voxels, voxels_by_roi,tumor_voxels, arcs, nodes, neighborhood, beamlets, deposition, proportionality,a,b):

    options = {
        "WLSACCESSID": "9a6728f5-f620-4829-b3b7-5e01b8e667d9",
        "WLSSECRET": "6f73eab7-cea6-4f6f-a167-510fd4aa9aa9",
        "LICENSEID": 2544834,
    }
    #env = Env(params=options)
    model = Model("optimal")#,env=env)  # create model

    m = 1000

    x = model.addVars(beamlets, vtype=GRB.CONTINUOUS, name="x")
    y = model.addVars(nodes, vtype=GRB.BINARY, name="y")
    z = model.addVars(voxels, vtype=GRB.BINARY, name="z")
    w = model.addVars([(node, voxel) for node in nodes for voxel in voxels], vtype= GRB.CONTINUOUS, name="w")
    s_pos = model.addVars(voxels, vtype=GRB.CONTINUOUS,name="s_pos")
    s_neg = model.addVars(voxels, vtype=GRB.CONTINUOUS,name="s_neg")
    #d[v] = sum(deposition[(b,v)]*x[b] for b in beamlets)
    peak_dose = 10
    valley_dose = 0


    model.addConstrs(y[i] + y[j] <= 1 for i,j in arcs)

    model.addConstrs(y[i] + sum(y[j] for j in neighborhood[i]) >= 1 for i in nodes)

    model.addConstrs(sum(deposition[(b, v)] * x[b] for b in beamlets) - ((peak_dose - valley_dose)*sum(proportionality[n][v] for n in nodes) + valley_dose)*z[v] - valley_dose*(1-z[v]) + s_pos[v] - s_neg[v] == 0 for v in voxels )

    model.addConstrs(z[v] >= sum(y[u] for u in {key: voxels[key] for key in neighborhood[v] if key in neighborhood[v]}.keys()) for v in nodes)

    for (n, v) in [(n, v) for n in nodes for v in voxels]:
        model.addConstr(w[n, v] <= z[v], name=f"Constraint_5f_{n}_{v}")

    # Equation (5g): w_{nv} <= p_{nv} for all n in N, v in V
    for (n, v) in [(n, v) for n in nodes for v in voxels]:
        model.addConstr(w[n, v] <= proportionality[n][v], name=f"Constraint_5g_{n}_{v}")

    # Equation (5h): w_{nv} >= p_{nv} - M * (1 - z_v) for all n in N, v in V
    for (n, v) in [(n, v) for n in nodes for v in voxels]:
        model.addConstr(w[n, v] >= proportionality[n][v] - m * (1 - z[v]), name=f"Constraint_5h_{n}_{v}")

    model.setObjective(quicksum(s_neg[v] + s_pos[v] for v in voxels.keys()), GRB.MINIMIZE)
    model.setParam("OutputFlag", 1)
    model.Params.TimeLimit = 21600
    #model.Params.TimeLimit = 160
    model.update()
    model.optimize()
    if model.status == GRB.INFEASIBLE:
        # Compute the IIS (Irreducible Infeasible Subsystem)
        model.computeIIS()
        # Output the IIS to a file
        model.write("infeasible_model.ilp")
        print("IIS written to infeasible_model.ilp")
    model.write("model.mps")

    # Optionally, increase the level of detail in the log
    model.setParam('OutputFlag', 1)  # Enable detailed output during optimization
    model.setParam('LogFile', 'gurobi_log.txt')
    x_values = {}
    y_values = {}
    p_values = {}
    z_t_values = {}
    z_o_values = {}
    z_neg_values = {}
    # Retrieve optimized values and split them into respective dictionaries
    for v in model.getVars():
        if v.VarName.startswith("x["):
            x_values[v.VarName] = v.X  # Use getAttr to get the value
        elif v.VarName.startswith("y["):
            y_values[v.VarName] = v.X  # Use getAttr to get the value
        elif v.VarName.startswith("p["):
            p_values[v.VarName] = v.X  # Store value in p_values
        elif v.VarName.startswith("z["):
            z_t_values[v.VarName] = v.X  # Use getAttr to get the value
        elif v.VarName.startswith("s_pos["):
            z_o_values[v.VarName] = v.X  # Store value in p_values
        elif v.VarName.startswith("s_neg["):
            z_neg_values[v.VarName] = v.X  # Store value in p_values
    # Return the dictionaries for each variable type
    return x_values, y_values, p_values,z_t_values,z_o_values, z_neg_values


voxels_by_roi = {}
with open("../fluence_map_models/data/voxels_by_roi.pkl", 'rb') as openfile:
    voxels_by_roi = pickle.load(openfile)

voxels = {}
with open("../fluence_map_models/data/voxels.pkl", 'rb') as openfile:
    voxels = pickle.load(openfile)

arcs = []
with open("../fluence_map_models/data/arcs.pkl", 'rb') as openfile:
    arcs = pickle.load(openfile)

nodes = {}
with open("../fluence_map_models/data/nodes.pkl", 'rb') as openfile:
    nodes = pickle.load(openfile)

neighborhood = {}
with open("../fluence_map_models/data/neighborhood.pkl", 'rb') as openfile:
    neighborhood = pickle.load(openfile)

beamlets = {}
with open("../fluence_map_models/data/beamlets.pkl", 'rb') as openfile:
    beamlets = pickle.load(openfile)

deposition = {}
with open("../fluence_map_models/data/deposition.pkl", 'rb') as openfile:
    deposition = pickle.load(openfile)

proportionality = {}
with open("../fluence_map_models/data/proportionality.pkl", 'rb') as openfile:
    proportionality = pickle.load(openfile)

tumor_voxels = {key: voxels[key] for key in voxels_by_roi["PTV GRID"] if key in voxels}

intensities, centers, penalties,target_deposition,z_o, s_neg = optimization2(voxels, voxels_by_roi, tumor_voxels, arcs, nodes, neighborhood, beamlets, deposition, proportionality, [-1,1],[10,-10])

optimal_nodes = []
'''
 for i in range(len(centers)):
    centers[i] = int(centers[i])
    if centers[i] == 1.0:
        optimal_nodes.append(centers.get(str(i)))
'''
print(intensities)
print({key: intensities[key] for key in intensities if intensities[key] > 0})
print(centers)
print({key: centers[key] for key in centers if centers[key] > .5})
print(len({key: centers[key] for key in centers if centers[key] > .5}))
print(penalties)
print({key: penalties[key] for key in penalties if penalties[key] < .5})
print(target_deposition)
print(z_o)
print(s_neg)
