import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# Read the CSV file (assumes no header)
df = pd.read_csv('dose_values.csv', header=None)

# Extract the values (second column) and reshape into a 10x10 grid
#values = df[1].values.reshape((10, 10))
values = df[1].values.reshape((10, 10))
# Plot the heatmap
plt.figure(figsize=(6, 6))
plt.imshow(values, cmap='hot', interpolation='nearest')
plt.colorbar(label='Gy')
plt.show()


df = pd.read_csv('deviation_values.csv', header=None)

# Extract the values (second column) and reshape into a 10x10 grid
values = df[1].values.reshape((10, 10))

# Plot the heatmap
plt.figure(figsize=(6, 6))
plt.imshow(values, cmap='hot', interpolation='nearest')
plt.colorbar(label='Gy')
plt.show()


# Define the voxel dictionary

voxel_dict = {'funny': [90, 80, 70, 60, 50, 40, 30], 'light': [1, 2, 3, 11, 12, 13, 21, 22, 23, 31, 32, 33, 41, 42, 43, 51, 52, 53, 61, 62, 63], 'plum': [91, 81, 71, 92, 93, 82, 83, 72, 73], 'tumor': [37, 38, 45, 46, 47, 48, 54, 55, 56, 57, 58, 64, 65, 66, 67, 68, 74, 75, 76, 77, 78, 85,86], 'border': [27,28, 49, 59,69,79, 35, 36, 39, 73, 44, 83, 84, 53, 95, 96, 87, 89, 88, 63]}


# Create a 10x10 grid initialized with zeros
grid = np.zeros((10, 10), dtype=int)

# Map each category to a unique integer
label_map = {
    "light": 1,   # gold
    "plum": 2,    # orchid
    "tumor": 3,   # deep red
    "funny": 4,   # slateblue
    "border": 5   # light red
}

# Fill the grid with category labels
for label, indices in voxel_dict.items():
    for idx in indices:
        row = idx // 10
        col = idx % 10
        grid[row, col] = label_map[label]

# Correct color alignment: index matches label_map value
colors = [
    'white',      # 0 - background
    'gold',       # 1 - light
    'orchid',     # 2 - plum
    '#cc0000',    # 3 - tumor (deep red)
    'slateblue',  # 4 - funny
    '#ff6666'     # 5 - border (light red)
]
cmap = ListedColormap(colors)

# Plot the heatmap
plt.figure(figsize=(6, 6))
im = plt.imshow(grid, cmap=cmap, interpolation='nearest')
plt.clim(0, 5)  # Adjust clim to match the label range (0 to 5)

# Add the colorbar with custom labels
'''
cbar = plt.colorbar(im, ticks=[0, 1, 2, 3, 4, 5])
cbar.set_label('Structure')
cbar.set_ticklabels([
    'No Structure', 'Organ 1 ', 'Organ 2',
    'Tumor', 'Organ 3', 'Tumor Border'
])
'''
#plt.title('Structure Map')
plt.grid(False)
plt.xticks(np.arange(10))
plt.yticks(np.arange(10))
plt.show()

exit()

# Flat list of labels (replace this with reading from a file if needed)
labels = [
    3,3,3,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
]

# Convert to NumPy array and reshape to 10x10
grid = np.array(labels).reshape((10, 10))

# Define a custom colormap for 0-3
colors = ['white', 'lightblue', 'violet', 'red']  # 0 to 3
cmap = ListedColormap(colors)

# Plot the heatmap
plt.figure(figsize=(6, 6))
plt.imshow(grid, cmap=cmap, interpolation='nearest')
plt.colorbar(ticks=[0, 1, 2, 3], label='Class Label')
plt.clim(0, 3)
plt.title('Class Label Heatmap')
plt.xlabel('Column')
plt.ylabel('Row')
plt.show()

# Full dataset pasted directly (tab-delimited)
data_str = """
3	1	0	0	0	0	0	0	0	0
3	1	0	0	0	0	0	0	0	0
3	1	0	0	0	0	0	0	0	0
2	1	0	0	0	0	0	0	0	0
2	1	0	0	0	0	0	0	0	0
2	1	0	0	0	0	0	0	0	0
1	1	0	0	0	0	0	0	0	0
1	1	0	0	0	0	0	0	0	0
1	1	0	0	0	0	0	0	0	0
1	1	0	0	0	0	0	0	0	0
1	3	1	0	0	0	0	0	0	0
1	3	1	0	0	0	0	0	0	0
1	3	1	0	0	0	0	0	0	0
1	2	1	0	0	0	0	0	0	0
1	2	1	0	0	0	0	0	0	0
1	2	1	0	0	0	0	0	0	0
1	1	1	0	0	0	0	0	0	0
1	1	1	0	0	0	0	0	0	0
1	1	1	0	0	0	0	0	0	0
1	1	1	0	0	0	0	0	0	0
0	1	3	1	0	0	0	0	0	0
0	1	3	1	0	0	0	0	0	0
0	1	3	1	0	0	0	0	0	0
0	1	2	1	0	0	0	0	0	0
0	1	2	1	0	0	0	0	0	0
0	1	2	1	0	0	0	0	0	0
0	1	1	1	0	0	0	0	0	0
0	1	1	1	0	0	0	0	0	0
0	1	1	1	0	0	0	0	0	0
0	1	1	1	0	0	0	0	0	0
0	0	1	3	1	0	0	0	0	0
0	0	1	3	1	0	0	0	0	0
0	0	1	3	1	0	0	0	0	0
0	0	1	2	1	0	0	0	0	0
0	0	1	2	1	0	0	0	0	0
0	0	1	2	1	0	0	0	0	0
0	0	1	1	1	0	0	0	0	0
0	0	1	1	1	0	0	0	0	0
0	0	1	1	1	0	0	0	0	0
0	0	1	1	1	0	0	0	0	0
0	0	0	1	3	1	0	0	0	0
0	0	0	1	3	1	0	0	0	0
0	0	0	1	3	1	0	0	0	0
0	0	0	1	2	1	0	0	0	0
0	0	0	1	2	1	0	0	0	0
0	0	0	1	2	1	0	0	0	0
0	0	0	1	1	1	0	0	0	0
0	0	0	1	1	1	0	0	0	0
0	0	0	1	1	1	0	0	0	0
0	0	0	1	1	1	0	0	0	0
0	0	0	0	1	3	1	0	0	0
0	0	0	0	1	3	1	0	0	0
0	0	0	0	1	3	1	0	0	0
0	0	0	0	1	2	1	0	0	0
0	0	0	0	1	2	1	0	0	0
0	0	0	0	1	2	1	0	0	0
0	0	0	0	1	1	1	0	0	0
0	0	0	0	1	1	1	0	0	0
0	0	0	0	1	1	1	0	0	0
0	0	0	0	1	1	1	0	0	0
0	0	0	0	0	1	3	1	0	0
0	0	0	0	0	1	3	1	0	0
0	0	0	0	0	1	3	1	0	0
0	0	0	0	0	1	2	1	0	0
0	0	0	0	0	1	2	1	0	0
0	0	0	0	0	1	2	1	0	0
0	0	0	0	0	1	1	1	0	0
0	0	0	0	0	1	1	1	0	0
0	0	0	0	0	1	1	1	0	0
0	0	0	0	0	1	1	1	0	0
0	0	0	0	0	0	1	3	1	0
0	0	0	0	0	0	1	3	1	0
0	0	0	0	0	0	1	3	1	0
0	0	0	0	0	0	1	2	1	0
0	0	0	0	0	0	1	2	1	0
0	0	0	0	0	0	1	2	1	0
0	0	0	0	0	0	1	1	1	0
0	0	0	0	0	0	1	1	1	0
0	0	0	0	0	0	1	1	1	0
0	0	0	0	0	0	1	1	1	0
0	0	0	0	0	0	0	1	3	1
0	0	0	0	0	0	0	1	3	1
0	0	0	0	0	0	0	1	3	1
0	0	0	0	0	0	0	1	2	1
0	0	0	0	0	0	0	1	2	1
0	0	0	0	0	0	0	1	2	1
0	0	0	0	0	0	0	1	1	1
0	0	0	0	0	0	0	1	1	1
0	0	0	0	0	0	0	1	1	1
0	0	0	0	0	0	0	1	1	1
0	0	0	0	0	0	0	0	1	3
0	0	0	0	0	0	0	0	1	3
0	0	0	0	0	0	0	0	1	3
0	0	0	0	0	0	0	0	1	2
0	0	0	0	0	0	0	0	1	2
0	0	0	0	0	0	0	0	1	2
0	0	0	0	0	0	0	0	1	1
0	0	0	0	0	0	0	0	1	1
0	0	0	0	0	0	0	0	1	1
0	0	0	0	0	0	0	0	1	1
"""


# Convert string to 2D NumPy array
lines = data_str.strip().splitlines()
data = np.array([[int(val) for val in line.split()] for line in lines])

# Reshape the data to create a 10x10 grid per column (10 columns, each a 10x10 image)
reshaped_data = data.reshape((10, 10, 10))  # 10 rows, 10 columns per dataset

# Use only red color maps for all beamlets
colormaps = ['Reds'] * 10

# Plot the data
fig, axes = plt.subplots(2, 5, figsize=(20, 10))  # 2 rows x 5 columns of subplots

for col in range(data.shape[1]):
    ax = axes[col // 5, col % 5]
    grid_data = reshaped_data[:, :, col]

    # Show the image using the red colormap
    cax = ax.imshow(grid_data, cmap='Reds', interpolation='nearest')
    ax.set_title(f"Beamlet {col + 1}", fontsize=12)
    fig.colorbar(cax, ax=ax)

plt.tight_layout()
plt.show()