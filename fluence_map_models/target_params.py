import numpy as np
import pandas as pd
import pickle
import time
from datetime import datetime

def is_point_in_sphere(point, center, radius=7.5):
    return np.linalg.norm(point - center) <= radius

def voxel_sample_proportion(center, voxel_pos, voxel_dims=(1, 1, 0.25), samples_per_dim=3):
    """
    Estimate the fraction of the voxel's volume that lies within a sphere centered at 'center'.
    """
    x_size, y_size, z_size = voxel_dims
    num_points_in_sphere = 0
    total_points = samples_per_dim ** 3

    # Get the lower corner of the voxel
    origin = np.array(voxel_pos) - np.array([x_size, y_size, z_size]) / 2

    # Generate sample points inside the voxel
    for i in range(samples_per_dim):
        for j in range(samples_per_dim):
            for k in range(samples_per_dim):
                offset = np.array([
                    (i + 0.5) * x_size / samples_per_dim,
                    (j + 0.5) * y_size / samples_per_dim,
                    (k + 0.5) * z_size / samples_per_dim,
                ])
                sample_point = origin + offset
                if is_point_in_sphere(sample_point, center):
                    num_points_in_sphere += 1

    return num_points_in_sphere / total_points

def get_target_dose_matrix_nonparallel(nodes, voxels, tumor_voxels, all_voxels, target_dose=15):
    matrix = []
    count = 0
    numtum = 0
    voxel_dims = (5,5,3)#(1, 1, 0.25)  # (x, y, z)
    for node, center in nodes.items():
        holder = []
        center = np.array(center)
        for voxel in range(all_voxels[-1][-1]):
            if voxel in voxels:
                numtum += 1
                voxel_pos = np.array(voxels[voxel]["position_xyz_mm"])
                proportion = voxel_sample_proportion(center, voxel_pos, voxel_dims=voxel_dims)
                dose = target_dose * proportion
                #holder.append(dose)
                holder.append(proportion)
                if proportion > 0:
                    count += 1
            else:
                holder.append(0)
        print("node: " + str(node))
        print("tumor voxels found: " + str(numtum))
        numtum = 0
        print("voxels with non-zero dose: " + str(count))
        count = 0
        matrix.append(holder)

    print("Calculated target dose matrix with shape " + str(len(matrix)) + ", " + str(len(matrix[0])))
    return matrix





if __name__ == "__main__":


    start_time = time.time()
    voxels_by_roi = {}
    with open("../fluence_map_models/data/voxels_by_roi.pkl", 'rb') as openfile:
        voxels_by_roi = pickle.load(openfile)

    nodes = {}
    with open("../fluence_map_models/data/nodes.pkl", 'rb') as openfile:
        nodes = pickle.load(openfile)

    voxels = {}
    with open("../fluence_map_models/data/voxels.pkl", 'rb') as openfile:
        voxels = pickle.load(openfile)

    data = []
    with open("../fluence_map_models/data/inf_params.pkl", 'rb') as openfile:
        data = pickle.load(openfile)

    voxel_coordinates, voxels_data, tumor_voxels, all_voxels, structure_list = data

    target_matrix = get_target_dose_matrix_nonparallel(nodes, voxels_data, all_voxels[structure_list.index("PTV")], all_voxels)

    with open('../fluence_map_models/data/target_matrix.pkl', 'wb') as outfile:
        pickle.dump(target_matrix, outfile)

    print("completed in " + str(time.time() - start_time))