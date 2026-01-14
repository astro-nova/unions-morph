
# 🌌 unions-morph

### Galaxy Cluster Dynamics & Evolution

[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Jupyter](https://img.shields.io/badge/Jupyter-Notebooks-orange.svg)](https://jupyter.org/)

*A project to study the dynamics and evolution of galaxy clusters using UNIONS data and simulations for the Taylor group*

[Overview](#overview) • [Structure](#project-structure) • [Workflow](#workflow) • [Contact](#contact)

---

</div>

## 🔭 Overview

This project investigates the connection between **dynamical history**, **galaxy orbits**, and their **structural properties** in galaxy clusters. By combining observational data from UNIONS (Ultraviolet Near Infrared Optical Northern Survey) with simulations, we aim to understand how galaxies evolve as they interact with the cluster environment.

## ✅ Work Done So Far

- ✨ **Sample Selection**: Galaxy samples for morphological analysis selected in `nbs/01 Morphology Example`
- 📐 **Morphology Analysis**: Comprehensive morphological parameters calculated and stored in `catalogs/morph.csv`

📓 Example Notebooks

| 📁 Notebook | Description |
|----------|-------------|
| `nbs/s01 Morphology Example.ipynb` | Sample selection, tile data fetching, and morphology calculation for all galaxies in a tile |


## Project Structure

```
unions-morph/
├── README.md
├── catalogs/           # Data catalogs
│   ├── morph.csv       # Morphology measurements
│   └── processed_tiles.csv # which tiles are analyzed (temporary)
├── lib/               # Core library modules
│   ├── io.py         # Input/output utilities
│   └── plotting.py   # Visualization tools
├── nbs/              # Jupyter notebooks
└── scripts/          # Processing scripts
    └── morph_parallel.py
```

**1️⃣ Clone the repository:**

```bash
git clone https://github.com/yourusername/unions-morph.git
cd unions-morph
```

**2️⃣ Create your own branch:**
```bash
git checkout -b your-branch-name
```

**3️⃣ Make changes and commit:**
```bash
git add .
git commit -m "Describe your changes"
git push origin your-branch-name
```

** 4️⃣ Open a pull request** in GitHub to merge your changes into the main branch.it push origin your-branch-name

When ready, open a pull request in github to merge your changes into the main branch.

## 📚 Citation

If you use this code or data in your research, please cite the Taylor group's relevant publications.

## 💬 Contact

For questions or collaboration inquiries, please contact the Taylor group.

