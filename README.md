# Shape Analysis via Skeletal Models

This repository contains software for fitting, analyzing, and visualizing the skeletal model **Locally Parameterized Discrete Skeletal Representations (LPDSRep)** of 3D anatomical objects.

LPDSRep is an invariant skeletal shape representation designed for statistical shape analysis. Unlike conventional discrete skeletal representations (DSRep), LPDSRep employs hierarchical local coordinate systems to represent skeletal geometry independently of global rigid transformations. This allows reliable statistical analysis without the bias introduced by alignment procedures.

The methodology was introduced in:

> Mohsen Taheri and Jörn Schulz.
> **Statistical Analysis of Locally Parameterized Shapes.**
> *Journal of Computational and Graphical Statistics*, 32(2):658–670, 2023.

---

## Synthetic Data Generation

This repository uses the **ETRep** R package for both model fitting and statistical analysis. In addition to analyzing real anatomical datasets, ETRep provides tools for generating synthetic tubular objects represented as **Elliptical Tube Representations (ETReps)**. The synthetic data are generated within the ETRep shape space using intrinsic transformations that preserve the geometric validity of the objects by enforcing the **Relative Curvature Condition (RCC)**. This enables the creation of realistic shape populations with controlled geometric variability while avoiding local self-intersections. 

The synthetic datasets are primarily intended for validating LP-DS-Rep fitting algorithms, evaluating statistical methods, and benchmarking shape analysis pipelines under controlled conditions. The generated ETReps can be converted into meshes, point clouds, and skeletal representations, making them suitable for simulation studies, hypothesis testing, machine learning, and visualization. The implementation follows the methodology introduced in the publication below.

> Mohsen Taheri Shalmani, Stephen M. Pizer, and Jörn Schulz.
> **The Mean Shape under the Relative Curvature Condition.**
> *Journal of Computational and Graphical Statistics*, 2025.
> https://doi.org/10.1080/10618600.2025.2535600

## Features

- LPDSRep model fitting
- Construction of hierarchical local coordinate systems
- Transformation between DSRep and LPDSRep
- Shape deformation and simulation
- Mean shape estimation
- Statistical hypothesis testing
- Shape visualization
- Mesh generation
- Shape correspondence analysis
- R implementation of the LPDSRep framework

---

## Workflow

The typical workflow consists of

1. Import an s-rep.
2. Convert the s-rep to LPDSRep.
3. Perform statistical analysis.
4. Compute mean shapes.
5. Simulate new shapes.
6. Visualize the results.

---

## Main Concepts

LPDSRep represents each skeletal atom using intrinsic geometric quantities defined relative to a local coordinate system.

The representation separates

- skeletal topology
- local orientations
- spoke directions
- spoke lengths
- local deformations

which enables interpretable analysis of local shape changes such as

- bending
- twisting
- stretching
- local expansion
- local contraction

without requiring global alignment.

---

## Applications

The framework has been developed primarily for medical image analysis and statistical shape analysis, including

- slab-shaped anatomical structures such as the hippocampus
- skeletal representations
- longitudinal shape analysis
- disease characterization

---

## Dependencies

The examples use

- R (≥ 4.2)
- ETRep
- rgl
- Rvcg
- Morpho
- shapes
- plotly
- ggplot2

---

## Citation

If you use this repository in your research, please cite

```bibtex
@article{Taheri2023LPDSRep,
  author = {Mohsen Taheri and J\"orn Schulz},
  title = {Statistical Analysis of Locally Parameterized Shapes},
  journal = {Journal of Computational and Graphical Statistics},
  volume = {32},
  number = {2},
  pages = {658--670},
  year = {2023},
  doi = {10.1080/10618600.2022.2116445}
}
```

---

## License

This project is released under the MIT License.

---

## Authors

**Mohsen Taheri Shalmani**
