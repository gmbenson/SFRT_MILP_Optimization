import math
import pickle
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import cupy as cp
from concurrent.futures import ProcessPoolExecutor
import matplotlib.pyplot as plt
import plotly.express as px

def get_tumor_bounds(voxels):
    ptv_voxels = voxels
    bounds = []
    for i in range(len(ptv_voxels[0])):
        bounds.append([max(list(map(list, zip(*ptv_voxels)))[i]), min(list(map(list, zip(*ptv_voxels)))[i])])

    return bounds


def gram_schmidt(vectors):
    """
    Orthonormalize a list of vectors using the Gram-Schmidt process.

    Args:
    vectors: List of linearly independent vectors.

    Returns:
    orthonormal_basis: List of orthonormal vectors.
    """
    # Create an empty list to hold the orthonormal basis
    orthonormal_basis = []

    for v in vectors:
        # Start with the original vector
        v = np.array(v, dtype=np.float64)

        # Subtract the projection of v onto each of the orthonormal basis vectors already found
        for b in orthonormal_basis:
            v -= np.dot(v, b) * b

        # Normalize the vector (make it a unit vector)
        v /= np.linalg.norm(v)

        # Append the orthonormal vector to the list
        orthonormal_basis.append(v)

    return np.array(orthonormal_basis)
def get_beamlets(voxels, bounds, beam = [], distance = 100, beamlet_count = 0):
    '''
    :param voxels:
    :param bounds:
    :param beam: list of [x,y,z] vector in direction of beam
    :return:
    '''

    beamlets = {}
    center = np.mean(bounds, axis = 1)
    beam = np.array(beam)
    beam = beam / np.linalg.norm(beam)
    arbitrary_vector = np.array([0,1,0])

    vector1 = np.cross(beam, arbitrary_vector)

    if np.linalg.norm(vector1) < 1e-6:
        arbitrary_vector = np.array([0, 0, 1])
        vector1 = np.cross(beam, arbitrary_vector)

    vector1 = vector1 / np.linalg.norm(vector1)
    vector2 = np.cross(beam, vector1)
    vector2 = vector2 / np.linalg.norm(vector2)

    #orthonormal_basis = gram_schmidt([vector1,vector2])
    num_beamlets = 5
    beamlet_spacing = 1
    count = beamlet_count
    for i in range(-1*num_beamlets, num_beamlets+1):
        for j in range(-1*num_beamlets, num_beamlets+1):
            offset = i *vector1*beamlet_spacing + j*vector2*beamlet_spacing
            beamlets[count] = center + offset + distance * beam
            count +=1

    return beamlets





def get_dose_deposition(tumor_voxels, voxels_by_roi, voxels, beams):
    """
    :param tumor_voxels: List of tumor voxel positions
    :param voxels_by_roi: Dictionary of voxels segmented by regions of interest
    :param voxels: Dictionary of all voxel positions
    :param beams: List of [x, y, z] unit vectors for each beam direction
    :return: beamlets and deposition dictionary
    """

    def fluence(d, r, D_0=1000, d_peak=10.0, sigma=30.0):
        return D_0 / ((d - d_peak) ** 2 + sigma) * np.exp(-(r / sigma) ** 2)

    beamlets = {}
    deposition = {}
    tumor_bounds = get_tumor_bounds(tumor_voxels)

    # Precompute constants
    beam_angle_norm_cache = {}

    # Cache beam direction norms outside the loops
    def get_beam_norm(beam):
        if tuple(beam) not in beam_angle_norm_cache:
            beam_angle_norm_cache[tuple(beam)] = np.linalg.norm(beam)
        return beam_angle_norm_cache[tuple(beam)]

    # Helper function to process each beam
    def process_beam(beam):
        # Normalize the beam direction
        beam = np.array(beam) / get_beam_norm(beam)  # Ensure it's a unit vector
        beamlet_holder = get_beamlets(voxels=tumor_voxels, bounds=tumor_bounds, beam=beam)
        print(f"Got beamlets for beam: {beam}")

        beamlet_results = {}

        for b, beamlet_origin in beamlet_holder.items():
            for v in voxels:
                voxel_loc_relative_to_beamlet = np.array(v) - beamlet_origin

                # Calculate distances and fluence
                rad = np.linalg.norm(np.cross(voxel_loc_relative_to_beamlet, beam)) / get_beam_norm(beam)
                dist = np.dot(voxel_loc_relative_to_beamlet, beam) / get_beam_norm(beam)
                beamlet_results[(int(b), int(v))] = fluence(dist, rad)

        return beamlet_holder, beamlet_results

    # Use ThreadPoolExecutor to parallelize beam processing
    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(process_beam, beam) for beam in beams]

        # Collect the results of each future
        for future in futures:
            beamlet_holder, beamlet_deposition = future.result()
            beamlets.update(beamlet_holder)
            deposition.update(beamlet_deposition)
            print(f"Got deposition for beam")

    return beamlets, deposition

def compute_fluence(v, beamlet, beam_angle):
    voxel_loc_relative_to_beamlet = cp.array(v) - beamlet
    # Calculate distance from beamlet
    rad = cp.linalg.norm(cp.cross(voxel_loc_relative_to_beamlet, beam_angle)) / cp.linalg.norm(beam_angle)
    dist = cp.dot(voxel_loc_relative_to_beamlet, beam_angle) / cp.linalg.norm(beam_angle)
    fluence = lambda d, r, D_0=1000, d_peak=10.0, sigma=30.0: D_0 / ((d - d_peak)**2 + sigma) * cp.exp(-(r / sigma)**2)
    return (int(beamlet), int(v), fluence(d=dist, r=rad))

def get_dose_deposition1(tumor_voxels, voxels_by_roi, voxels, beams):
    """
    :param voxels:
    :param beams: list including angles of each beam
    :return:
    """

    beamlets = {}
    deposition = {}
    beamlet_number = 0
    tumor_bounds = get_tumor_bounds(tumor_voxels)
    beamlet_count = 0

    for beam in beams:
        beamlet_holder = get_beamlets(voxels=tumor_voxels, bounds= tumor_bounds, beam = beam, beamlet_count = beamlet_count)
        print("got beamlets for beam:", beam)
        for b in beamlet_holder:
            for v in voxels:
                voxel_loc_relative_to_beamlet = np.array(v) - beamlet_holder[b]
                #calculate distance from beamlet
                beam_angle = np.array(beam)

                rad = np.linalg.norm(np.cross(voxel_loc_relative_to_beamlet, beam_angle)) / np.linalg.norm(beam_angle)
                dist = np.dot(voxel_loc_relative_to_beamlet, beam_angle) / np.linalg.norm(beam_angle)
                fluence = lambda d, r, D_0=1000, d_peak=10.0, sigma=30.0: D_0 / ((d - d_peak)**2 + sigma) * np.exp(-(r / sigma)**2)
                deposition[(int(b),int(v))] = fluence(d = dist,r = rad)
        beamlets.update(beamlet_holder)
        print("got deposition for beam:", beam)
    return beamlets, deposition



voxels = {}
with open("../fluence_map_models/data/voxels.pkl", 'rb') as openfile:
    voxels = pickle.load(openfile)

voxels_by_roi = {}
with open("../fluence_map_models/data/voxels_by_roi.pkl", 'rb') as openfile:
    voxels_by_roi = pickle.load(openfile)

nodes = {}
with open("../fluence_map_models/data/nodes.pkl", 'rb') as openfile:
    nodes = pickle.load(openfile)

tumor_voxels = [voxels[key] for key in voxels_by_roi["PTV GRID"] if key in voxels]

beamlets, deposition = get_dose_deposition(voxels= voxels,voxels_by_roi = voxels_by_roi, tumor_voxels=tumor_voxels, beams = [[1,0,0],[0,0,1]])

with open('../fluence_map_models/data/beamlets.pkl', 'wb') as outfile:
    pickle.dump(beamlets, outfile)

with open('../fluence_map_models/data/deposition.pkl', 'wb') as outfile:
    pickle.dump(deposition, outfile)

print("done")