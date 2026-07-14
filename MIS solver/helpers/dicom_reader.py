import pydicom
import json

def extract_candidate_points_from_dicom(
    dicom_path,
    roi_name='GTV-1',
    granularity=10,
    output_json_path='../data/test_candidate_points.json',
    method=None
):
    if method is None:
        raise ValueError("A discretization method (callable) must be provided")

    ds = pydicom.dcmread(dicom_path)
    structures = {item.ROINumber: item.ROIName for item in ds.StructureSetROISequence}

    try:
        k = list(structures.values()).index(roi_name)
    except ValueError:
        raise ValueError(f"ROI '{roi_name}' not found in StructureSetROISequence")

    contours = ds.ROIContourSequence

    x, y, z = [], [], []

    def get_planar_points(slice_idx):
        contour_seq = contours[k].ContourSequence[slice_idx]
        z_val = contour_seq.ContourData[2]

        points = []
        i = 0
        while i < len(contour_seq.ContourData):
            points.append([contour_seq.ContourData[i], contour_seq.ContourData[i + 1]])
            i += 3

        # Use the passed method (discretization function) here:
        points_holder = method(points=points, granularity=granularity, plot=False)

        candidates = []
        for p in points_holder:
            candidates.append([p[0], p[1], float(z_val)])
            x.append(p[0])
            y.append(p[1])
            z.append(z_val)
        return candidates

    candidate_points = []
    for i in range(len(contours[k].ContourSequence)):
        candidate_points.append(get_planar_points(i))
        #print(f'Processed slice: {i+1}/{len(contours[k].ContourSequence)}')

    #print(f"Finished placing {len(x)} points across {len(contours[k].ContourSequence)} slices")

    # Save to JSON
    with open(output_json_path, 'w') as outfile:
        json.dump(candidate_points, outfile)

    return candidate_points



if __name__ == "__main__":
    dicom_file = r"C:\Users\Grant\Desktop\manifest-1739568518674\NSCLC-Radiomics\LUNG1-006\01-01-2014-StudyID-NA-99263\4.000000-NA-27681\1-1.dcm"
    roi = "GTV-1"
    gran = 10
    output_path = '../data/test_candidate_points005.json'
    discretization_module = 'grid_placement'  # or 'relative_grid_placement', etc.

    extract_candidate_points_from_dicom(
        dicom_path=dicom_file,
        roi_name=roi,
        granularity=gran,
        discretization_module_name=discretization_module,
        output_json_path=output_path
    )
