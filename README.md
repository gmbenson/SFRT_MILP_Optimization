# Mathematical Optimization for Spatially Fractionated Radiation Therapy (SFRT)

> Master's thesis project combining graph theory, mixed-integer optimization, and radiation treatment planning to optimize Spatially Fractionated Radiation Therapy (SFRT).

---

## Overview

Spatially Fractionated Radiation Therapy (SFRT) delivers radiation through a pattern of high-dose "peaks" separated by lower-dose "valleys", allowing aggressive treatment of large tumors while reducing damage to healthy tissue.

Current treatment planning largely relies on manual sphere placement by clinicians.

This project formulates the entire planning process as a mathematical optimization problem.

The optimization pipeline simultaneously determines

- where radiation spheres should be placed,
- how many spheres should be used,
- and the optimal beamlet intensities required to deliver the prescribed dose while satisfying clinical dose constraints.

The resulting framework combines discrete optimization, convex optimization, computational geometry, and medical image processing into a single optimization workflow.

---

# Optimization Pipeline

```text
CT Scan
   │
   ▼
Tumor Contours (ROI)
   │
   ▼
Candidate Sphere Generation
   │
   ▼
Graph Construction
   │
   ▼
Maximum Independent Set
   │
   ▼
Dose Prescription Matrix
   │
   ▼
Mixed Integer Fluence Optimization
   │
   ▼
Deliverable Radiation Treatment Plan
```

---

# Major Components

## 1. Tumor Discretization

Designed a computational geometry algorithm that converts physician-defined tumor contours into a discrete set of feasible radiation sphere locations.

Features include

- point-in-polygon detection
- ray casting
- boundary interpolation
- tangent handling
- support for disconnected tumor regions
- configurable grid discretization
- square and triangular sampling strategies

---

## 2. Graph Construction

Each feasible sphere center becomes a graph node.

Edges connect sphere centers that violate the minimum separation distance required by clinical SFRT guidelines.

Finding the optimal sphere placement becomes a Maximum Independent Set problem.

Implemented and compared three formulations:

- Arc-based MILP
- Node-based MILP
- Quadratic Binary Optimization

---

## 3. Fluence Map Optimization

Extended classical fluence map optimization by coupling beamlet intensity optimization with sphere placement.

Instead of optimizing a fixed prescription matrix, the prescription itself becomes a function of the optimization variables.

The optimization simultaneously determines

- sphere locations
- beamlet intensities
- delivered voxel dose

while satisfying

- voxel dose limits
- mean dose constraints
- DVH constraints
- organ-at-risk constraints

---

## 4. Large Scale Optimization

Patient treatment plans contain

- hundreds of thousands of voxels
- tens of thousands of beamlets
- multi-gigabyte sparse dose matrices

Several techniques were required to keep optimization tractable:

- CSR sparse matrices
- voxel downsampling
- warm-started branch-and-bound
- L1 objective linearization
- preprocessing pipelines
- proportional voxel overlap approximation

---

# Engineering Challenges

| Challenge | Problem | Solution |
|---|---|---|
| Irregular 3D tumor geometry | Radiation sphere placement requires identifying all feasible locations inside complex, patient-specific tumor contours. | Developed a custom CT-based discretization pipeline using point-in-polygon methods, ray casting, boundary interpolation, and geometric filtering to generate candidate sphere locations. |
| Combinatorial explosion of sphere placement | The number of possible sphere configurations grows exponentially as tumor size and discretization resolution increase. | Formulated sphere placement as a Maximum Independent Set (MIS) problem on a conflict graph and solved using mixed-integer optimization techniques. |
| Large-scale graph optimization | Patient-specific tumor discretization can generate thousands of candidate sphere locations and millions of potential conflicts. | Implemented graph construction algorithms and compared multiple MIS formulations to identify efficient optimization strategies. |
| Coupling discrete decisions with radiation delivery | Optimal sphere placement alone does not guarantee a clinically deliverable treatment plan. | Integrated MIS-based sphere selection with fluence map optimization to jointly optimize sphere locations and beamlet intensities. |
| Multi-gigabyte dose matrices | Radiation dose deposition matrices can exceed 10 GB due to the number of voxels and beamlets required for patient-specific planning. | Utilized sparse matrix representations (CSR format), preprocessing, and voxel downsampling to make large-scale optimization computationally feasible. |
| High-dimensional mixed-integer optimization | Simultaneously optimizing sphere placement, beamlet intensities, and clinical constraints creates a challenging optimization problem. | Applied MILP formulations, L1 objective linearization, branch-and-bound optimization, and warm-start strategies to improve solver performance. |
| Accurate dose prescription modeling | Traditional fluence optimization assumes a fixed target dose distribution, which does not capture variable SFRT sphere placement. | Developed a dynamic dose prescription matrix where the target dose changes based on optimization-selected sphere locations. |
| Computational geometry approximation | Exact sphere-voxel intersection calculations are computationally expensive for patient-scale models. | Developed a proportional voxel overlap approximation to estimate dose contribution while maintaining tractable computation times. |
| Clinical constraint integration | Treatment plans must satisfy safety constraints for healthy organs while maintaining tumor dose requirements. | Incorporated voxel limits, mean dose constraints, and dose-volume histogram (DVH) constraints into the optimization model. |

---

# Technologies

### Optimization

- Mixed Integer Programming (MILP)
- Linear Programming
- Convex Optimization
- Maximum Independent Set
- Branch-and-Bound
- Multi-objective Optimization

### Scientific Computing

- Python
- NumPy
- SciPy
- Sparse Linear Algebra

### Medical Computing

- PortPy
- CT Scan Processing
- DICOM-derived ROI data
- Dose-volume histogram (DVH) constraints

### Solvers

- Gurobi

---

# Research Contributions

This project introduces several novel components beyond standard radiation treatment planning.

### Graph-based sphere placement

Developed a graph representation of feasible SFRT sphere locations generated directly from CT-derived tumor contours.

### Coupled optimization formulation

Integrated Maximum Independent Set optimization with Fluence Map Optimization, allowing sphere placement and beamlet intensities to be optimized simultaneously.

### Dynamic dose prescription

Introduced a variable dose prescription matrix whose values depend on optimization decisions rather than remaining fixed throughout optimization.

### Practical optimization improvements

Implemented

- sparse matrix preprocessing
- warm starts
- objective linearization
- voxel overlap approximation

to solve clinically-sized optimization problems.

---

# Repository Structure

```
src/
    preprocessing/
    graph_generation/
    optimization/
    fluence/
    visualization/

data/
    patient_data/

results/
    treatment_plans/

figures/
```

---

# Example Workflow

1. Load CT scan using PortPy
2. Extract tumor contours
3. Generate candidate sphere centers
4. Construct incompatibility graph
5. Solve Maximum Independent Set
6. Build dose prescription matrix
7. Compute dose deposition matrix
8. Solve mixed-integer fluence optimization
9. Export treatment plan

---

# Future Work

Potential extensions include

- VMAT optimization
- Robust optimization under patient motion
- GPU acceleration
- Multi-criteria optimization
- Adaptive radiation therapy

Description goes here...

This repository supports analysis of dicom structural files

* Note: all path variables must be changed to the relevant CT\RS\ file

## Installation

To run, download poetry with `pip install poetry`.

Once poetry is downloaded, configure virtual environment with `poetry install`.

After virtual environment has been created, activate with `poetry shell`.

### Running Ingest

To run ingest, make sure you are in the root directory `Research/`, and then execute `python ingest.py`.

## Description

### primary pipeline

dicom_reader.py - creates a list of all points within a given ROI, prints to json file - requires grid_placement.py - writes candidate_points.json (slices of node lists)

parameter_writer.py - reads candidate_points.json, calculates nodes and arcs, writes node_data.json, arc_data.json, arc_list.json, adjacency_matrix.json files

Graph_Optimization.py - a preliminary implementation of the Maximum Independent Set problem in Gurobi, reads node_data.json and arc_data.json

**All optimization files require dicom_reader.py and parameter_writer.py

## Secondary Models

#### Other Optimization Models:

* Clique_Optimization.py : solves MIS problem using all maximal cliques rather than edges. Needs cliques.json from Clique_Problem.py
    
* Node_Optimization.py : solves MIS problem using aggregated constraints, such that there is only one adjacency constraint per node. Needs arc_list.json and node_data.json
To initialize env, run `poetry install`

#### Heuristics:

* cube_heuristic.py : finds a solution for the MIS problem by assuming that the spheres all sit on a uniform grid, which is optimal for a rectangular shape. Finds every possible grid solution and returns best solution.
    
* Greedy_Algorithm.py : Implementation of a Greedy algorithm that performs suspiciously poorly. Needs validation, potentially correction.
 
    

## Other Files

* Clique_Problem.py - uses networkX library to find all cliques and create cliques.json. Requires node_data.json and arc_data.json

* grid_placement.py - containes helper funcitons for dicom_reader.py

* Structure_Visualizer.py - can be used to visualize ROIs in 3D and 2D slices



## Miscellaneous files

* cube_optimization.py - another implementation of the MIS problem to test computation time, data storage
* dicom_tester.py - reads dicom files, don't use
