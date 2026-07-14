'''

Given a RStruct DICOM file, extract contours and plot an indexed ROI as a 3d plot and a single slice of the ROI

'''
import os
import pydicom
import csv
import matplotlib.pyplot as plt
import plotly.express as px
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

df = []

directory = r"C:\Users\Grant\PycharmProjects\SFRT_Optimization\misc\portpy_dataset"
for file in os.listdir("portpy_dataset"):
    f = os.path.join(directory, file)
    ds = pydicom.dcmread(f)
    if ds.Modality != 'RTSTRUCT':
        print(f"File {f} is not an RTSTRUCT file.")

    contours = ds.ROIContourSequence

    structures = {}
    for item in ds.StructureSetROISequence:
       structures[item.ROINumber] = item.ROIName

    points = []

    k = list(structures.values()).index('GTV-1')

    #create lists for plotting values
    x = []
    y = []
    z = []

    for j in range(len(contours[k].ContourSequence)):
       i = 0
       while i < len(contours[k].ContourSequence[j].ContourData):
           #add (x,y,z) contour point to list points
           points.append([contours[k].ContourSequence[j].ContourData[i], contours[k].ContourSequence[j].ContourData[i+1], contours[k].ContourSequence[j].ContourData[i+2]])
           x.append(contours[k].ContourSequence[j].ContourData[i])
           y.append(contours[k].ContourSequence[j].ContourData[i+1])
           z.append(contours[k].ContourSequence[j].ContourData[i+2])
           i+=3
    volume = 0
    width = [max(x) - min(x), max(y) - min(y), max(z) - min(z)]
    result = [ds.PatientID]
    result.extend(width)
    df.append(result)
    print(df)

filename = "PortPy_Dataset_Tumor_Size.csv"
with open(filename, mode = "w", newline='') as file:
    writer = csv.writer(file)
    for row in df:
        writer.writerow(row)

print(f'CSV file "{filename}" has been created successfully.')


