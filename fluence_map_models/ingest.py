import subprocess
from voxel_ingest import run_voxel_extraction
from parameter_ingest import run_fluence_pipeline

# Define your input/output
rtstruct_path = r"C:\Users\Grant\Desktop\manifest-1739568518674\NSCLC-Radiomics\LUNG1-046\01-01-2014-StudyID-NA-84290\2.000000-NA-06055\1-1.dcm"
output_dir = r"C:\Users\Grant\Desktop\output_voxels"
rois = ["GTV-1"]

print("running node extraction")
# Run the function
run_voxel_extraction(rtstruct_path, output_dir, rois)

# You can change this to any patient ID from your PortPy dataset
patient_id = "Lung_Patient_27"
data_dir = r"C:\Users\Grant\Downloads\PortPy_data"

print("running fluence calculation")

run_fluence_pipeline(patient_id=patient_id, data_dir=data_dir)


print("running target parameter calculation")
subprocess.run(["python", "target_params.py"], check=True)
subprocess.run(["python", "matrix_sparsification.py"], check= True)