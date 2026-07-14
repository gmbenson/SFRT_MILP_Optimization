'''

Given a RStruct DICOM file, extract contours and plot an indexed ROI as a 3d plot and a single slice of the ROI

'''

import pydicom
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import plotly.express as px
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


#load file r"path"
ds = pydicom.dcmread(r"C:\Users\Grant\Desktop\manifest-1739568518674\NSCLC-Radiomics\LUNG1-042\01-01-2014-StudyID-NA-81603\4.000000-NA-19421\1-1.dcm")
#ds = pydicom.dcmread(r"C:\Users\Grant\Downloads\Anonymized grid\Grid 2 anonymized\2024-07__Studies\Grid2_Grid2_RTst_2024-07-29_145646_._ARIA.RadOnc.Structure.Sets_n1__00000\2.16.840.1.114362.1.12306304.27066498827.682660245.301.301.dcm")

patient = 34
#extract contour names
contours = ds.ROIContourSequence

structures = {}
for item in ds.StructureSetROISequence:
   structures[item.ROINumber] = item.ROIName

print(structures.values())


points = []
#k = 12#39  #ROI 12 = ptv_grid 19 = ptv spheres, 23 = PTV VMAT, 27 GRIDptv_TM_RESEARCH
i = 0  #data point
j = 2  #slice



k = list(structures.values()).index('GTV-1') ##'Grid Boundary 0.8 cm Inner Margin')
#k = list(structures.values()).index('Grid_Boundary_SR')

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


# Create the 3D plot
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot points
scatter = ax.scatter3D(x, y, z, c=z, cmap='viridis', marker='o', s=20, edgecolor='k', alpha=0.8)

# Add colorbar (optional, especially useful if 'z' has meaningful variation)
#cbar = fig.colorbar(scatter, ax=ax, pad=0.1)
#cbar.set_label('Z Value')

# Set axis labels with formatting
ax.set_xlabel('X (mm)', fontsize=12, labelpad=10)
ax.set_ylabel('Y (mm)', fontsize=12, labelpad=10)
ax.set_zlabel('Z (mm)', fontsize=12, labelpad=10)

# Set title
ax.set_title(f'Patient {patient}', fontsize=14, pad=20)

# Improve tick label size
ax.tick_params(axis='both', which='major', labelsize=10)

# Set view angle (optional, can adjust for better perspective)
ax.view_init(elev=20, azim=135)

# Tight layout
plt.tight_layout()
plt.show()
exit()

ifig = px.scatter_3d(x = x, y = y, z=z)
ifig.show()

print(str(len(points)) + " border/contour points in the ROI")

x1 = []
y1 = []

j=15 #slice
i=0
while i < len(contours[k].ContourSequence[j].ContourData):
    points.append([contours[k].ContourSequence[j].ContourData[i], contours[k].ContourSequence[j].ContourData[i + 1]])
    x1.append(contours[k].ContourSequence[j].ContourData[i])
    y1.append(contours[k].ContourSequence[j].ContourData[i + 1])
    i += 3

fig2 = plt.figure()
ax = plt.axes()
ax.scatter(x1,y1)
plt.show()