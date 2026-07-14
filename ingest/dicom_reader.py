import pydicom
import json
from grid_placement import get_candidate_points

def extract_candidate_points_from_dicom(
    dicom_path,
    roi_name='GTV-1',
    granularity=10,
    output_json_path='../data/test_candidate_points.json'
):
    """
    Given an RStruct DICOM file path and ROI name,
    extract internal candidate points for each contour slice.

    Returns:
      candidate_points: list of slices, each slice a list of [x, y, z] points
    Also saves JSON to output_json_path.
    """
    ds = pydicom.dcmread(dicom_path)
    structures = {item.ROINumber: item.ROIName for item in ds.StructureSetROISequence}

    # Find index k corresponding to ROI name
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

        points_holder = get_candidate_points(points=points, granularity=granularity, plot=False)
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
        print(f'Processed slice: {i+1}/{len(contours[k].ContourSequence)}')

    # Optional: validate shapes
    for slice_pts in candidate_points:
        for pt in slice_pts:
            if len(pt) != 3:
                print("Warning: point with incorrect dimension found:", pt)

    print(f"Finished placing {len(x)} points across {len(contours[k].ContourSequence)} slices")

    # Save to JSON
    with open(output_json_path, 'w') as outfile:
        json.dump(candidate_points, outfile)

    return candidate_points


if __name__ == "__main__":
    # Example usage - replace these paths as needed
    dicom_file = r"C:\Users\Grant\Desktop\manifest-1739568518674\NSCLC-Radiomics\LUNG1-006\01-01-2014-StudyID-NA-99263\4.000000-NA-27681\1-1.dcm"
    roi = "GTV-1"
    gran = 10
    output_path = '../data/test_candidate_points005.json'

    extract_candidate_points_from_dicom(
        dicom_path=dicom_file,
        roi_name=roi,
        granularity=gran,
        output_json_path=output_path
    )
