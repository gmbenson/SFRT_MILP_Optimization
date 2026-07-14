import pydicom
import pickle
from fluence_map_models.grid_placement import get_candidate_points
import itertools
import os
import SimpleITK as sitk
from rt_utils import RTStructBuilder
from scipy.ndimage import binary_erosion
from rt_utils import RTStructBuilder
import time
from datetime import datetime
import scipy.ndimage as ndimage
import numpy as np

def get_planar_points(j, roi_index, contours, granularity = 1):
    '''
    Takes contours, read from dicom file, outputs a list of x,y,z pairs
    :param j:
    :param contours:
    :param granularity:
    :return:
    '''

    i = 0
    z_val = contours[roi_index].ContourSequence[j].ContourData[2]
    points = []
    candidates = []

    while i < len(contours[roi_index].ContourSequence[j].ContourData):
        points.append([contours[roi_index].ContourSequence[j].ContourData[i], contours[roi_index].ContourSequence[j].ContourData[i + 1]])
        i += 3
    points_holder = get_candidate_points(points = points,granularity= granularity, plot = False)
    for i in range(len(points_holder)):
        candidates.append([points_holder[i][0], points_holder[i][1], float(z_val)])
    return candidates



def extract_voxels_from_structure(contours, roi_index):

    points_holder = []
    for i in range(len(contours[roi_index].ContourSequence)):
        points_holder.append(get_planar_points(j=i, roi_index=roi_index, contours=contours, granularity=2.5))  # granularity in (mm)
        print('slice: ' + str(i))

    #deslice
    points_list = []
    for i in range(len(points_holder)):
        for j in range(len(points_holder[i])):
            points_list.append(points_holder[i][j])
    return points_list
def extract_nodes_from_structure(contours, roi_index):

    points_holder = []
    for i in range(len(contours[roi_index].ContourSequence)):
        points_holder.append(get_planar_points(j=i, roi_index=roi_index, contours=contours, granularity=10))  # granularity in (mm)
        print('slice: ' + str(i))

    #deslice
    points_list = []
    for i in range(len(points_holder)):
        for j in range(len(points_holder[i])):
            points_list.append(points_holder[i][j])
    return points_list


def get_voxels(path = None, rois = []):
    if path == None:
        ds = pydicom.dcmread(r"C:\Users\Grant\Downloads\RS.1.2.246.352.205.4769015767137605396.12843964432630046891.dcm")
    else:
        ds = pydicom.dcmread(path)

    contours = ds.ROIContourSequence

    structures = {}
    for item in ds.StructureSetROISequence:
        structures[item.ROINumber] = item.ROIName
    print(structures.values())

    voxels = {}
    voxels_by_roi = {}
    voxel_num = 0
    '''
    for roi in rois:
        roi_index = list(structures.values()).index(roi)
        voxel_holder = extract_voxels_from_structure(contours=contours, roi_index = roi_index)
        starting_voxel_num = voxel_num
        for i in voxel_holder:
            voxels[voxel_num] = i
            voxel_num += 1
        voxels_by_roi[roi] = [i for i in range(starting_voxel_num,voxel_num-1)]
        print("extracted voxels for " + roi)
    '''
    nodes = {}
    tumor_index = list(structures.values()).index("GTV-1")
    node_holder = extract_nodes_from_structure(contours=contours, roi_index=tumor_index)
    node_num = 0
    for i in node_holder:
        nodes[node_num] = i
        node_num += 1
    print("extracted", node_num, "nodes")
    '''
    deposition = {}
    x1 = [0,-2.5]
    y1 = [0,-2.5]
    z1 = [0,-2.5]
    combin1 = list(itertools.product(x1,y1,z1))

    combin2 = [[0,0,5]]
    for index, coord in nodes.items():
        if index not in deposition:
            deposition[index] = {}
        for offset in combin1:
            for key, value in voxels.items():
                if value[0] == coord[0] + offset[0] and value[1] == coord[1] + offset[1] and value[2] == coord[1] + offset[2]:
                    deposition[index][key] = 1
        for offset in combin2:
            deposition[index] = {key:0.75 for key,value in voxels.items() if value[0] == coord[0] + offset[0] and value[1] == coord[1] + offset[1] and value[2] == coord[2] + offset[2]}
    x1 = [2.5,]


    print(deposition)
    '''
    return voxels, voxels_by_roi, nodes

def get_spacing_from_dicom(ct_dir):
    for file in os.listdir(ct_dir):
        if file.endswith(".dcm"):
            ds = pydicom.dcmread(os.path.join(ct_dir, file))
            spacing_xy = ds.PixelSpacing
            spacing_z = ds.SliceThickness
            return [float(spacing_xy[0]), float(spacing_xy[1]), float(spacing_z)]
    raise RuntimeError("No valid DICOM CT slices found in the directory.")

def preprocess(ct_dir, rtstruct_path, structure_name="GTV-1", output_path="rtstruct_shrunk.dcm"):
    # Load RTSTRUCT in context of CT series
    builder = RTStructBuilder.create_from(
        dicom_series_path=ct_dir,
        rt_struct_path=rtstruct_path
    )

    # Get mask for structure (e.g. GTV-1)
    mask = builder.get_roi_mask_by_name(structure_name)

    # Get voxel spacing from CT slice metadata
    spacing = get_spacing_from_dicom(ct_dir)

    # Shrink ROI by 8 mm in all directions
    shrink_mm = 8
    radius_vox = [int(shrink_mm / s) for s in spacing]
    radius = max(radius_vox)  # Conservative shrink across all axes

    # Perform binary erosion
    eroded_mask = binary_erosion(mask, iterations=radius)

    # Add new ROI and save modified RTSTRUCT
    builder.add_roi(mask=eroded_mask.astype(bool), name=f"PTV_shrunk_{shrink_mm}mm")

    builder.save(output_path)


#preprocess(r"C:\Users\Grant\Desktop\manifest-1739568518674\NSCLC-Radiomics\LUNG1-005\01-01-2014-StudyID-NA-93819\1.000000-NA-68747", r"C:\Users\Grant\Desktop\manifest-1739568518674\NSCLC-Radiomics\LUNG1-005\01-01-2014-StudyID-NA-93819\2.000000-NA-99068\1-1.dcm")

rois = ["PTV_shrunk_8mm"]#, "Brain"]#, "Brain", "BrachialPlexs", "Esophagus","Eye_L","Eye_R","Larynx","SpinalCord","Teeth","Thyroid"]

def run_voxel_extraction(rtstruct_path, output_dir, rois=["PTV_shrunk_8mm"]):
    import os
    os.makedirs(output_dir, exist_ok=True)
    print("extracting nodes")

    voxels, voxels_by_roi, nodes = get_voxels(rois=rois, path=rtstruct_path)

    with open(f"data/nodes.pkl", 'wb') as outfile:
        pickle.dump(nodes, outfile)

    print("nodes loaded for ", rois)


    maxb = [-1e10, -1e10, -1e10]
    minb = [1e10, 1e10, 1e10]
    for coord in nodes.values():
        for i in range(3):
            minb[i] = min(minb[i], coord[i])
            maxb[i] = max(maxb[i], coord[i])
    print("node min bounds: " + str(minb))
    print("node max bounds: " + str(maxb))

if __name__ == "__main__":
    rtstruct_path = r"C:\Users\Grant\Desktop\manifest-1739568518674\NSCLC-Radiomics\LUNG1-016\01-01-2014-StudyID-NA-36302\4.000000-NA-72678\1-1.dcm"
    output_dir = r"C:\Users\Grant\Desktop\output_voxels"
    rois = ["GTV-1"]

    print("running voxel extraction")
    # Run the function
    run_voxel_extraction(rtstruct_path, output_dir, rois)