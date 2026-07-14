import pickle
import math
import random
import portpy.photon as pp
import gc
import subprocess
import time
from datetime import datetime

def create_arcs(nodes, min_dist = 30):
    arcs = []
    node_arcs = {}
    matrix = []
    neighborhood = {key: [] for key in nodes}
    num_arcs = 0
    for i in nodes.keys():
        count = 0

        holder = []
        holder.clear()
        matrix_holder = []
        matrix_holder.clear()
        for j in nodes.keys():
            if math.dist(nodes[i], nodes[j]) <= min_dist and i != j:
                arcs.append([i,j])
                neighborhood[i].append(j)
                num_arcs +=1
    return arcs, node_arcs, matrix, neighborhood, num_arcs


def node_voxel_intersect_test(node_center, radius, voxel_center, voxel_size_xy, voxel_size_z):
    """
    Check if the bounding box of a voxel intersects the sphere.
    """
    min_dist = math.sqrt(2*(voxel_size_xy**2) + voxel_size_z**2)
    if radius - (0.5*min_dist) <= math.dist(node_center, voxel_center) <= radius + (0.5*min_dist):
        return 0
    elif radius - (0.5*min_dist) >= math.dist(node_center, voxel_center):
        return -1
    elif radius + (0.5*min_dist) <= math.dist(node_center, voxel_center):
        return 1
    else:
        raise Exception("NOOOOOOOOOOO")

def exact_proportion(node, radius, voxel_center,voxel_size_xy,voxel_size_z, num_samples = 1000):
    x_c, y_c, z_c = node
    x_v, y_v, z_v = voxel_center


    inside_count = 0
    for i in range(num_samples):
        # Generate a random point within the voxel
        x_rand = random.uniform(x_v - voxel_size_xy*0.5, x_v + voxel_size_xy*0.5)
        y_rand = random.uniform(y_v - voxel_size_xy*0.5, y_v + voxel_size_xy*0.5)
        z_rand = random.uniform(z_v - voxel_size_z*0.5, z_v + voxel_size_z*0.5)

        # Calculate the distance from the point to the center of the sphere
        distance = math.sqrt((x_rand - x_c) ** 2 + (y_rand - y_c) ** 2 + (z_rand - z_c) ** 2)

        # If the point is inside the sphere, increase the count
        if distance <= radius:
            inside_count += 1

    # Calculate the proportion of the voxel inside the sphere
    proportion_inside = inside_count / num_samples
    return proportion_inside



def node_deposition_to_voxels(key, node, radius, voxels, voxel_size_xy, voxel_size_z):
    sphere_proportions = {}
    for voxel in voxels:
        intersection = node_voxel_intersect_test(node, radius, voxels[voxel],voxel_size_xy,voxel_size_z)
        if intersection == 0:
            proportion = exact_proportion(node, radius, voxels[voxel],voxel_size_xy,voxel_size_z)
            sphere_proportions[voxel] = proportion
        elif intersection == 1:
            sphere_proportions[voxel] = 0
        elif intersection == -1:
            sphere_proportions[voxel] = 1
        else:
            raise Exception("Gosh darn it!")
    return (key, sphere_proportions)



def get_portpy_data(data_dir =r"C:\Users\Grant\Downloads\PortPy_data", patient_id = "Lung_Patient_4"):


    #data_dir = r"C:\Users\Grant\Downloads\Lung_Patient_15"
    data = pp.DataExplorer(data_dir=data_dir)
    data.patient_id = patient_id

    ct = pp.CT(data)
    stucts = pp.Structures(data)

    beams = pp.Beams(data)
    clinical_criteria = pp.ClinicalCriteria(data, protocol_name="Lung_2Gy_30Fx")
    opt_vox_res = [ct_res * factor for ct_res, factor in zip(ct.get_ct_res_xyz_mm(), [5,5,1])]
    inf_matrix = pp.InfluenceMatrix(ct=ct, structs=stucts, beams=beams)
    inf_matrix_dv = inf_matrix.create_down_sample(opt_vox_xyz_res_mm=opt_vox_res)
    plan = pp.Plan(structs=stucts,beams = beams, inf_matrix=inf_matrix_dv)

    voxel_coordinates = {}
    A_matrix = inf_matrix_dv.A.toarray()

    print("finished calculating fluence matrix")

    voxel_data = {}
    count = 0
    tumor_voxels = []

    #for j in range(len(plan.structures.opt_voxels_dict["voxel_idx"])):
    #for i in plan.structures.opt_voxels_dict["voxel_idx"][plan.structures.opt_voxels_dict['name'].index('PTV')]:
    for i in inf_matrix_dv.opt_voxels_dict["voxel_idx"][inf_matrix_dv.opt_voxels_dict["name"].index("PTV")]:
        voxel_data[i] = plan.inf_matrix.get_voxel_info(i)
        voxel_coordinates[i] = voxel_data[i]["position_xyz_mm"]
        try:
            if len(voxel_data[i]) >= 2:
                if 'PTV' in voxel_data[i]["structures"]:
                    tumor_voxels.append(i)
                    #if count == 1000:
                    #    break
                    if count % 1000 == 0:
                        print(count)
            count += 1
        except KeyError:
            print("error")
    print("found " + str(count) + " PTV voxels")

    return A_matrix, voxel_coordinates, voxel_data, tumor_voxels, inf_matrix_dv.opt_voxels_dict["voxel_idx"], inf_matrix_dv.opt_voxels_dict["name"]


def run_fluence_pipeline(patient_id="Lung_Patient_4", data_dir=r"C:\Users\Grant\Downloads\PortPy_data"):
    start_time = time.time()

    # Load nodes
    with open("../fluence_map_models/data/nodes.pkl", 'rb') as openfile:
        nodes = pickle.load(openfile)

    inf_matrix, voxel_coordinates, voxels_data, tumor_voxels, all_voxels, structure_list = get_portpy_data(data_dir=data_dir, patient_id=patient_id)
    inf_params = [voxel_coordinates, voxels_data, tumor_voxels, all_voxels, structure_list]

    with open("../fluence_map_models/data/voxel_data", "wb") as outfile:
        pickle.dump(voxels_data, outfile)

    with open("../fluence_map_models/data/inf_params.pkl", 'wb') as outfile:
        pickle.dump(inf_params, outfile)

    arcs, node_arcs, matrix, neighborhood, num_arcs = create_arcs(nodes)

    print(len(arcs), "arcs loaded")

    with open('../fluence_map_models/data/arcs.pkl', 'wb') as outfile:
        pickle.dump(arcs, outfile)

    with open('../fluence_map_models/data/inf_matrix.pkl', 'wb') as outfile:
        pickle.dump(inf_matrix, outfile)

    print("completed in " + str(time.time() - start_time))

if __name__ == "__main__":
    run_fluence_pipeline(patient_id="Lung_Patient_11", data_dir=r"C:\Users\Grant\Downloads\PortPy_data")