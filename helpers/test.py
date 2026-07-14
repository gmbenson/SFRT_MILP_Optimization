import sys
import portpy.photon as pp
import os
import matplotlib.pyplot as plt
import cvxpy as cp
# specify the patient data location.
import pandas as pd
import h5py
import numpy as np

# specify the patient data location.
data_dir = r"C:\Users\Grant\Downloads\Lung_Patient_11"
# Use PortPy DataExplorer class to explore PortPy data
data = pp.DataExplorer(data_dir=data_dir)
# Pick a patient


f =  h5py.File(r"C:\Users\Grant\Downloads\Lung_Patient_11\Lung_Patient_11\CT_Data.h5", "r")

print(list(f.keys()))

dset = f["ct_hu_3d"]
print(dset.shape)


data.patient_id = 'Lung_Patient_11'

data.get_tcia_metadata()
exit()


ct = pp.CT(data)
structs = pp.Structures(data)

# Select beam ids based upon target location. We select beam ids from 0 to 36 based upon beams metadata
beam_ids = np.arange(0, 37) #  users can modify this based upon their beam selection
beams = pp.Beams(data, beam_ids=beam_ids)

# load clinical criteria from the config files for which plan to be optimized
protocol_name = 'Lung_2Gy_30Fx'
clinical_criteria = pp.ClinicalCriteria(data, protocol_name=protocol_name)

# Loading hyper-parameter values for optimization problem
protocol_name = 'Lung_2Gy_30Fx_vmat'
vmat_opt_params = data.load_config_opt_params(protocol_name=protocol_name)

# # Creating optimization structures (i.e., Rinds, PTV-GTV)
structs.create_opt_structures(opt_params=vmat_opt_params,
                              clinical_criteria=clinical_criteria)

# Loading influence matrix
inf_matrix = pp.InfluenceMatrix(ct=ct, structs=structs, beams=beams)

# Assign discrete beam/control_point_ids to arcs and create arcs dictionary.
# Below is an example of creating 2 arcs. Users can create single or multiple arcs.
arcs_dict = {'arcs': [{'arc_id': "01", "beam_ids": beam_ids[0]}]}
# Create an object of Arcs class
#arcs = pp.Arcs(arcs_dict=arcs_dict, inf_matrix=inf_matrix)

# Create a plan using ct, structures, beams and influence matrix. Clinical criteria is optional

my_plan = pp.Plan(ct=ct, structs=structs, beams=beams, inf_matrix=inf_matrix, clinical_criteria=clinical_criteria)#, arcs=arcs)
'''
vmat_opt = pp.VmatScpOptimization(my_plan=my_plan,
                                  opt_params=vmat_opt_params)
# Run Sequential convex algorithm for optimising the plan.
# The final result will be stored in sol and convergence will store the convergence history (i.e., results of each iteration)
sol, convergence = vmat_opt.run_sequential_cvx_algo(solver='MOSEK', verbose=True)

# Visualize convergence. The convergence dataframe contains the following columns:
df = pd.DataFrame(convergence, columns=['outer_iteration', 'inner_iteration', 'step_size_f_b', 'forward_backward', 'intermediate_obj_value', 'actual_obj_value', 'accept'])
# We can, for example, plot the actual and intermediate objective values against the outer iteration
df.plot(x='outer_iteration', y=['actual_obj_value', 'intermediate_obj_value'])
plt.show()

my_plan = pp.load_plan(plan_name='my_plan_vmat.pkl', path=os.path.join(r'C:\temp', data.patient_id))
sol = pp.load_optimal_sol(sol_name='sol_vmat.pkl', path=os.path.join(r'C:\temp', data.patient_id))
'''
out_rt_plan_file = r'rt_plan_portpy_vmat.dcm'  # change this file directory based upon your needs
in_rt_plan_file = r'rt_plan_echo_vmat.dcm'  # change this directory as per your
#pp.write_rt_plan_vmat(my_plan=my_plan, in_rt_plan_file=in_rt_plan_file, out_rt_plan_file=out_rt_plan_file)



# get the corresponding tcia collection/subject ID
print(data.get_tcia_metadata())
exit()
# Specify the location and name of the DICOM RT Dose file
dose_file_name = os.path.join(r'C:\temp', data.patient_id, 'rt_dose_portpy_vmat.dcm')
# Convert the DICOM dose into PortPy format
ecl_dose_3d = pp.convert_dose_rt_dicom_to_portpy(my_plan=my_plan, dose_file_name=dose_file_name)
ecl_dose_1d = inf_matrix.dose_3d_to_1d(dose_3d=ecl_dose_3d)

exit()











beams_df, structs_df = data.display_patient_metadata(return_beams_df=True, return_structs_df=True)
structs_df.head()

ct = pp.CT(data)
structs = pp.Structures(data)
beams = pp.Beams(data)

##### FLUENCE MATRIX ##########

inf_matrix = pp.InfluenceMatrix(ct=ct, structs=structs, beams=beams)


# Load ct, structure set, beams for the above patient using CT, Structures, and Beams classes
ct = pp.CT(data)
structs = pp.Structures(data)
beams = pp.Beams(data)
# Pick a protocol
protocol_name = 'Lung_2Gy_30Fx'
# Load clinical criteria for a specified protocol
clinical_criteria = pp.ClinicalCriteria(data, protocol_name=protocol_name)

# Load hyper-parameter values for optimization problem for a specified protocol
opt_params = data.load_config_opt_params(protocol_name=protocol_name)
# Create optimization structures (i.e., Rinds)
structs.create_opt_structures(opt_params=opt_params, clinical_criteria=clinical_criteria)
# Load influence matrix
inf_matrix = pp.InfluenceMatrix(ct=ct, structs=structs, beams=beams)

# Create a plan using ct, structures, beams and influence matrix, and clinical criteria
my_plan = pp.Plan(ct=ct, structs=structs, beams=beams, inf_matrix=inf_matrix, clinical_criteria=clinical_criteria)

print(my_plan.inf_matrix.get_voxel_info(0))
exit()

# Create cvxpy problem using the clinical criteria and optimization parameters
opt = pp.Optimization(my_plan, opt_params=opt_params, clinical_criteria=clinical_criteria)
opt.create_cvxpy_problem()
# Solve the cvxpy problem using Mosek
sol = opt.solve(solver='MOSEK', verbose=False)