import csv
import pickle
import numpy as np

matrix = []
'''

with open("../fluence_map_models/data/target_matrix.pkl", "rb") as openfile:
    matrix = pickle.load(openfile)

A = np.array(matrix)
B = [i.sum() for i in A]
'''

with open(r"C:\Users\Grant\Downloads\mis_constr_dose_results.csv", mode = 'r') as openfile:
    reader = csv.reader(openfile)
    count = 0
    full_count = 0
    for row in reader:
        full_count +=1
        if float(row[0]) > 0:
            print("node: " + str(count))
            count +=1
    roi = row[38968+2018+5930+691+18663+18461+21329:38968+2018+5930+691+18663+18461+21329+912]
    print(str(full_count)+ " nodes found")
    print("total nodes used: "+ str(count))