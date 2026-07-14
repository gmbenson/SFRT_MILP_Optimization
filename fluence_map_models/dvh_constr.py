import numpy as np
import pandas as pd
import pickle
import portpy as pp

def get_max_dose_constr(voxels_data, voxels, struct, dose):
    return 0
def get_dvh_constrs(nodes, voxels, all_voxels, target_dose = 15):
    matrix = []
    count = 0
    numtum = 0
    for node, center in nodes.items():
        holder = []
        center = np.array(center)
        for voxel in all_voxels[9]:
            if voxel in voxels.keys():
                numtum +=1
            #if 'PTV' in voxels[voxel]["structures"]:
                if np.linalg.norm(np.array(voxels[voxel]["position_xyz_mm"]) - center) <= 15: #(15^2)
                    holder.append(target_dose)
                    count += 1
                else:
                    holder.append(0)
            else:
                holder.append(0)
        print("node: "+ str(node))
        print("tumor voxels found: " + str(numtum))
        numtum = 0
        print("voxels inside nodes found: " + str(count))
        count = 0
        matrix.append(holder)

    print("found " + str(count) + " voxels within nodes")
    print("Calculated target dose matrix with shape " + str(len(matrix)) + ", " + str(len(matrix[0])))
    return matrix




if __name__ == "__main__":
    voxels_by_roi = {}
    with open("../fluence_map_models/data/voxels_by_roi.pkl", 'rb') as openfile:
        voxels_by_roi = pickle.load(openfile)

    nodes = {}
    with open("../fluence_map_models/data/nodes.pkl", 'rb') as openfile:
        nodes = pickle.load(openfile)

    voxels = {}
    with open("../fluence_map_models/data/voxels.pkl", 'rb') as openfile:
        voxels = pickle.load(openfile)

    data_dir = r"C:\Users\Grant\Downloads\Lung_Patient_15"
    data = pp.DataExplorer(data_dir=data_dir)
    data.patient_id = "Lung_Patient_15"
    clinical_criteria = pp.ClinicalCriteria(data, protocol_name="Lung_2Gy_30Fx")

    voxel_coordinates, voxels_data, tumor_voxels, all_voxels = data
    exit()
    target_matrix = get_dvh_constrs(nodes, voxels_data, all_voxels)

    with open('../fluence_map_models/data/target_matrix.pkl', 'wb') as outfile:
        pickle.dump(target_matrix, outfile)
