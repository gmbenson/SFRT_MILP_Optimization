import gc
import os
import importlib
from helpers.dicom_reader import extract_candidate_points_from_dicom
from helpers.parameter_writer import process_nodes
from models.arc_optimization import optimize_arc_mis
from models.node_optimization import optimize_node_mis
from models.adjacency_optimization import optimize_adjacency_mis
import sys

class Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for s in self.streams:
            s.write(data)
            s.flush()  # to ensure it shows immediately

    def flush(self):
        for s in self.streams:
            s.flush()

logfile = open("log_output.txt", "w", encoding="utf-8", buffering=1)
tee = Tee(sys.__stdout__, logfile)  # sys.__stdout__ is the original console
sys.stdout = tee
sys.stderr = tee  # optional: also log errors



# List of discretization methods to test
discretization_methods = [
    "grid_placement",
    "relative_grid_placement",
    "triangle_grid_placement",
    "relative_triangle_grid_placement",
]

DATA_FOLDER = "data"                  # Input: contains only valid DICOM files
OUTPUT_FOLDER = "results"            # Output: root directory for outputs
ROI_NAME = "PTV"
GRANULARITY = 10  # in mm

def get_discretization_func(module_name):
    module = importlib.import_module(f"algos.{module_name}")
    return getattr(module, "get_candidate_points")


def process_patient(dicom_fname):
    patient_id = os.path.splitext(dicom_fname)[0]
    dicom_path = os.path.join(DATA_FOLDER, dicom_fname)

    for method_name in discretization_methods:
        print(f"")
        print(f"")
        print(f"[{patient_id}] → {method_name}")
        print(f"")
        # Step 1: Extract candidate points
        discretizer = get_discretization_func(method_name)
        candidate_json = os.path.join(
            OUTPUT_FOLDER, patient_id, method_name, "candidate_points.json"
        )
        os.makedirs(os.path.dirname(candidate_json), exist_ok=True)

        candidate_points = extract_candidate_points_from_dicom(
            dicom_path=dicom_path,
            granularity=GRANULARITY,
            output_json_path=candidate_json,
            method=discretizer
        )

        candidate_points = [pt for slice_pts in candidate_points for pt in slice_pts]

        print(f"  ✓ {len(candidate_points)} candidate processed")

        # Step 2: Write node/arc/adjacency parameters
        params = process_nodes(candidate_points, OUTPUT_FOLDER)

        nodes = params["nodes"]
        arcs = params["arcs"]
        adjacency = params["adjacency"]
        num_nodes = params["num_nodes"]
        num_arcs = params["num_arcs"]
        neighborhood = params["neighborhood"]

        print(f"")
        print(f"  ✓ Graph built ({num_nodes} nodes, {num_arcs} edges)")
        print(f"")

        # Step 3: Arc MIS
        optimal_nodes_arc, mis_size_arc, gap, rt = optimize_arc_mis(
            nodes, arcs, time_limit=600, show_plot=False
        )
        print(f"  ✓ Arc MIS size: {mis_size_arc}, with gap {gap} for {rt} seconds")

        # Step 4: Node MIS
        optimal_nodes_node, mis_size_node, gap, rt = optimize_node_mis(
            nodes, neighborhood, time_limit=600, show_plot=False
        )
        print(f"  ✓ Node MIS size: {mis_size_node}, with gap {gap} for {rt} seconds")

        # Step 5: Adjacency MIS
        optimal_nodes_adj, mis_size_adj, gap, rt = optimize_adjacency_mis(
            nodes, adjacency, time_limit=600, show_plot=False
        )
        print(f"  ✓ Adjacency MIS size: {mis_size_adj}, with gap {gap} for {rt} seconds")

        gc.collect()  # Clear memory for next run

def main():
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    for fname in os.listdir(DATA_FOLDER):
        if fname.lower().endswith(".dcm"):
            print(f"")
            print(f"\n=== Processing {fname} ===")
            print(f"")
            process_patient(fname)

    print("\n✅ All done.")

if __name__ == "__main__":
    main()
