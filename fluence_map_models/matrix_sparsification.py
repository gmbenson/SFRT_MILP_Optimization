import pickle
import concurrent.futures
import os
import multiprocessing
from scipy.sparse import csr_matrix
import numpy as np
import os


# Define the process_row function at the module level (outside of create_list_np)
def process_row(i, matrix, progress_counter):
    holder = {}
    for j in range(len(matrix[i])):
        if matrix[i][j] != 0:
            holder[j] = matrix[i][j]
    # Directly update the progress_counter value without using get_lock
    progress_counter.value += 1
    return (i, holder)


def create_list(matrix):
    adj_list = {}

    # Use ThreadPoolExecutor or ProcessPoolExecutor, depending on your workload
    max_workers = os.cpu_count()  # You can adjust this number if needed

    with multiprocessing.Manager() as manager:
        progress_counter = manager.Value('i', 0)  # 'i' stands for integer

        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = []
            for i in range(len(matrix)):
                futures.append(executor.submit(process_row, i, matrix, progress_counter))

            # Process the results and print progress directly from the main thread
            for idx, future in enumerate(concurrent.futures.as_completed(futures)):
                i, holder = future.result()
                adj_list[i] = holder

                # Print progress every 100 rows
                if (idx + 1) % 100 == 0:
                    print(f"Processed {idx + 1}/{len(matrix)} rows.")

    return adj_list

def create_list_np(matrix):
    adj_list = {}
    for i in range(len(matrix)):
        holder = {}
        for j in range(len(matrix[i])):
            if matrix[i][j] != 0:
                holder[j] = matrix[i][j]
        adj_list[i] = holder
        if i % 1000 == 0:
            print(str(int(i/1000)) + '/' + str(round(len(matrix)/1000)))
    return adj_list

if __name__ == '__main__':
    data = []
    with open("../fluence_map_models/data/inf_params.pkl", 'rb') as openfile:
        data = pickle.load(openfile)

    voxel_pos, voxel_dat, tumor_voxel, all_voxels, structure_list = data
    matrix = []
    with open("../fluence_map_models/data/inf_matrix.pkl", 'rb') as openfile:
        matrix = pickle.load(openfile)

    target = []
    with open("../fluence_map_models/data/target_matrix.pkl", 'rb') as openfile:
        target = pickle.load(openfile)

    print("matrix loaded")
    '''
    adj_list = create_list_np(matrix)
    target_list = create_list_np(target)

    '''
    for row in target:
        row.extend([0] * (len(matrix) - len(row)))
    sparse_matrix = csr_matrix(matrix)
    target = np.array(target)
    target = np.transpose(target)
    target_matrix = csr_matrix(target)
    
    print("List created")

    with open("../fluence_map_models/data/sparse_matrix.pkl", 'wb') as openfile:
        pickle.dump(sparse_matrix, openfile)

    with open("../fluence_map_models/data/sparse_target_matrix.pkl", 'wb') as openfile:
        pickle.dump(target_matrix, openfile)

    #os.remove("../fluence_map_models/data/inf_matrix.pkl")
    #os.remove("../fluence_map_models/data/target_matrix.pkl")

    '''
    with open("../fluence_map_models/data/inf_list.pkl", 'wb') as openfile:
        pickle.dump(adj_list, openfile)
    with open("../fluence_map_models/data/target_list.pkl", 'wb') as openfile:
        pickle.dump(target_list, openfile)
    print("List saved")'''