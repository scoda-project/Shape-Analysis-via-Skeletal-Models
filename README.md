# LPDSRep: Local Parameterized Discrete Skeletal Representation

This repository contains software for fitting, analyzing, and visualizing **Locally Parameterized Discrete Skeletal Representations (LPDSRep)** of 3D anatomical objects.

LPDSRep is an invariant skeletal shape representation designed for statistical shape analysis. Unlike conventional discrete skeletal representations (DSRep), LPDSRep employs hierarchical local coordinate systems to represent skeletal geometry independently of global rigid transformations. This allows reliable statistical analysis without the bias introduced by alignment procedures.

The methodology was introduced in:

> Mohsen Taheri and Jörn Schulz.
> **Statistical Analysis of Locally Parameterized Shapes.**
> *Journal of Computational and Graphical Statistics*, 32(2):658–670, 2023.

---

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
