import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import json


# Function to read the CSV and plot a 3D scatterplot
def plot_3d_scatter(csv_file):
    # Read CSV file into a pandas DataFrame (without headers)
    data = pd.read_csv(csv_file, header=None)

    # Ensure there are exactly 3 columns
    if data.shape[1] != 3:
        raise ValueError("CSV must contain exactly 3 columns for x, y, and z coordinates.")

    # Extract x, y, z values
    x = data[0]  # First column
    y = data[1]  # Second column
    z = data[2]  # Third column

    # Create a figure and 3D axis
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Create a 3D scatter plot
    ax.scatter(x, y, z, c='r', marker='o')

    # Add labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Patient 8 Optimal Grid Placement')

    # Show the plot
    plt.show()


def load_json(filename):
    with open(filename, 'r') as file:
        return json.load(file)


# Function to plot 3D coordinates from the loaded JSON data
def plot_json_coordinates(file):
    data = load_json(file)
    # Create a 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Initialize empty lists to store x, y, and z values
    x_all = []
    y_all = []
    z_all = []

    # Loop through each list of coordinates
    for coordinates in data:
        # Assuming each entry in 'coordinates' is a list of [x, y, z] points
        x_all.extend([point[0] for point in coordinates])
        y_all.extend([point[1] for point in coordinates])
        z_all.extend([point[2] for point in coordinates])

    # Plot all points on the same scatter plot
    ax.scatter(x_all, y_all, z_all, c='b', marker='o')

    # Label the axes
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('3D Coordinates Plot')

    # Display the plot
    plt.show()

# Example usage
csv_file = r"C:\Users\Grant\PycharmProjects\SFRT_Optimization\data\candidate_points001.json"  # Replace with your CSV file path
plot_json_coordinates(csv_file)
