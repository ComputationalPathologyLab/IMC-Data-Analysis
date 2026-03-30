# IMC Data Analysis Pipeline

![IMC](https://img.shields.io/badge/Technology-IMC-blue)
![R](https://img.shields.io/badge/R-Bioconductor-green)
![Python](https://img.shields.io/badge/Python-3.x-yellow)
![Docker](https://img.shields.io/badge/Docker-Steinbock-blue)
![Status](https://img.shields.io/badge/Status-Active-success)

A comprehensive, reproducible workflow for **Imaging Mass Cytometry (IMC)** data analysis, covering preprocessing, segmentation, feature extraction, phenotyping, and spatial visualization using **Steinbock**, **Python**, and **Bioconductor (R)**.

---

# Table of Contents

- [1. Introduction](#1-introduction)
- [2. Imaging Mass Cytometry](#2-imaging-mass-cytometry)
- [3. Workflow Overview](#3-workflow-overview)
- [4. Input Data Structure](#4-input-data-structure)
- [5. Directory Preparation](#5-directory-preparation)
- [6. Channel Stacking](#6-channel-stacking)
- [7. Panel Annotation](#7-panel-annotation)
- [8. Steinbock Processing](#8-steinbock-processing)
- [9. Image Metadata](#9-image-metadata)
- [10. R Analysis](#10-r-analysis)
- [11. Marker Selection](#11-marker-selection)
- [12. Dimensionality Reduction and Clustering](#12-dimensionality-reduction-and-clustering)
- [13. Heatmap and Marker Interpretation](#13-heatmap-and-marker-interpretation)
- [14. Cell Type Annotation](#14-cell-type-annotation)
- [15. Spatial Visualization](#15-spatial-visualization)
- [16. Reproducibility](#16-reproducibility)
- [17. Repository Scope](#17-repository-scope)
- [18. Author](#18-author)

---

# 1. Introduction

This repository presents a complete and reproducible framework for the analysis of Imaging Mass Cytometry data. The workflow begins with channel-wise OME-TIFF files and proceeds through image stacking, segmentation, feature extraction, clustering, phenotyping, and spatial visualization.

The objective of the pipeline is to convert multiplexed tissue imaging data from the pixel level into biologically interpretable single-cell and spatial information. The framework is structured to support reproducibility, modularity, and adaptation to additional regions of interest (ROIs) and future datasets.

---

# 2. Imaging Mass Cytometry

Imaging Mass Cytometry is a high-dimensional spatial proteomics technology that enables simultaneous detection of multiple protein markers within a tissue section. Each marker is measured as an independent image channel, and pixel intensity reflects marker abundance at a given spatial location.

The raw data generally consist of multiple TIFF images per ROI, each representing one marker. At this stage, the data preserve spatial information but do not yet define cells, cell boundaries, cell identities, or cell-cell relationships. Downstream interpretation therefore requires systematic preprocessing and integration.

---

# 3. Workflow Overview

```text
Single-channel OME-TIFF files
        ↓
Channel stacking into one multi-channel TIFF
        ↓
Panel annotation
        ↓
Segmentation with Steinbock (DeepCell Mesmer)
        ↓
Cell-level intensity extraction
        ↓
Morphological feature extraction
        ↓
Spatial neighbor graph construction
        ↓
Export to CSV / H5AD / GraphML
        ↓
Import into R as SpatialExperiment
        ↓
Arcsinh transformation
        ↓
PCA, UMAP, clustering
        ↓
Cluster interpretation and cell type annotation
        ↓
Spatial visualization