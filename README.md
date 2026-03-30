# IMC Data Analysis Pipeline

![IMC](https://img.shields.io/badge/Technology-IMC-blue)
![R](https://img.shields.io/badge/R-Bioconductor-green)
![Python](https://img.shields.io/badge/Python-3.x-yellow)
![Docker](https://img.shields.io/badge/Docker-Steinbock-blue)
![Status](https://img.shields.io/badge/Status-Active-success)

A comprehensive, reproducible workflow for Imaging Mass Cytometry (IMC) data analysis, covering preprocessing, segmentation, feature extraction, phenotyping, and spatial visualization using Steinbock, Python, and Bioconductor (R).

---

# Table of Contents

1. Introduction  
2. Imaging Mass Cytometry  
3. Workflow Overview  
4. Input Data Structure  
5. Directory Preparation  
6. Channel Stacking  
7. Panel Annotation  
8. Steinbock Processing  
9. Image Metadata  
10. R Analysis  
11. Marker Selection  
12. Dimensionality Reduction  
13. Heatmap  
14. Cell Annotation  
15. Spatial Visualization  
16. Saving Outputs  
17. Troubleshooting  
18. Reproducibility  
19. Repository Scope  
20. Author  

---

# 1. Introduction

This repository presents a complete and reproducible framework for the analysis of Imaging Mass Cytometry data. The workflow begins with channel-wise OME-TIFF files and proceeds through image stacking, segmentation, feature extraction, clustering, phenotyping, and spatial visualization.

The objective of the pipeline is to convert multiplexed tissue imaging data from the pixel level into biologically interpretable single-cell and spatial information. The framework is structured to support reproducibility, modularity, and scalability across multiple regions of interest and datasets.

---

# 2. Imaging Mass Cytometry

Imaging Mass Cytometry is a high-dimensional spatial proteomics technology that enables simultaneous detection of multiple protein markers within a tissue section. Each marker is measured as an independent image channel, and pixel intensity reflects marker abundance at a given spatial location.

The raw data consist of multiple OME-TIFF files per ROI, each corresponding to a single marker channel. These images preserve spatial resolution but lack explicit cell-level segmentation, necessitating downstream computational processing.

---

# 3. Workflow Overview

```text
Single-channel OME-TIFF files
        ↓
Channel stacking into multi-channel TIFF
        ↓
Panel annotation
        ↓
Segmentation (DeepCell Mesmer)
        ↓
Cell-level intensity extraction
        ↓
Morphological feature extraction
        ↓
Spatial neighbor graph construction
        ↓
Export (CSV / H5AD / GraphML)
        ↓
Import into R (SpatialExperiment)
        ↓
Transformation and normalization
        ↓
Dimensionality reduction and clustering
        ↓
Cell type annotation
        ↓
Spatial visualization

This sequence reflects the logical transition from raw image data to biologically interpretable cellular and spatial outputs

# 4. Input Data Structure

A typical input region of interest (ROI) is organized as a directory containing multiple single-channel OME-TIFF files, where each file corresponds to one protein marker acquired during Imaging Mass Cytometry.

```text
ROI001_D13/
├── 115In_CD20.ome.tiff
├── 127I_127I.ome.tiff
├── 131Xe_131Xe.ome.tiff
├── 134Xe_134Xe.ome.tiff
├── 138Ba_138Ba.ome.tiff
├── 141Pr_CD68.ome.tiff
├── 142Nd_CD138.ome.tiff
├── 143Nd_CD47.ome.tiff
├── 144Nd_CD31.ome.tiff
├── ...
└── 193Ir_DNA2.ome.tiff
```
Each TIFF file represents a single acquisition channel corresponding to a specific metal-tagged antibody or control channel. The file naming convention typically encodes both the metal isotope and the biological marker.

At this stage, the dataset preserves high-resolution spatial information across all markers but lacks cell-level segmentation and quantitative summaries. Consequently, the data are not yet suitable for downstream single-cell or spatial analyses.

To enable compatibility with Steinbock and subsequent analytical frameworks, these individual channel images must be consolidated into a single multi-channel TIFF file while preserving channel order and spatial alignment.

---

# 5. Directory Preparation

Purpose: To establish a structured directory layout for preprocessing and downstream analysis.

The following command creates the directory required for storing multi-channel TIFF files:

```bash
mkdir -p data/img
```

Interpretation: The data/img/ directory serves as the image input location for Steinbock. Each ROI must be represented as a single multi-channel TIFF within this folder prior to segmentation and feature extraction. The standardized directory structure ensures compatibility with Steinbock commands and downstream processing pipelines.

# 6. Channel Stacking

Purpose: To convert single-channel TIFF files into a single multi-channel image compatible with Steinbock.

Steinbock requires each region of interest to be represented as a multi-channel TIFF file, where all marker channels are stacked in a consistent and predefined order. This transformation ensures that spatial alignment across channels is preserved and enables downstream segmentation and feature extraction.

The resulting image has the following structure:

```text
Channels × Height × Width
```
## 6.1 Channel Stacking Script

```python
from pathlib import Path
import sys
import tifffile as tiff
import numpy as np

# Usage:
# python stack_one_roi.py ROI001_D13 data/img/ROI001_D13.tiff

if len(sys.argv) >= 3:
    roi_dir = Path(sys.argv[1])
    out_path = Path(sys.argv[2])
else:
    roi_dir = Path("ROI001_D13")
    out_path = Path("data/img/ROI001_D13.tiff")

out_path.parent.mkdir(parents=True, exist_ok=True)

channel_order = [
    "115In_CD20.ome.tiff",
    "127I_127I.ome.tiff",
    "131Xe_131Xe.ome.tiff",
    "134Xe_134Xe.ome.tiff",
    "138Ba_138Ba.ome.tiff",
    "141Pr_CD68.ome.tiff",
    "142Nd_CD138.ome.tiff",
    "143Nd_CD47.ome.tiff",
    "144Nd_CD31.ome.tiff",
    "145Nd_Tbet.ome.tiff",
    "146Nd_BCL2.ome.tiff",
    "147Sm_CD44.ome.tiff",
    "148Nd_CD163.ome.tiff",
    "149Sm_CD45RO.ome.tiff",
    "150Nd_PDL1.ome.tiff",
    "151Eu_MHC_II.ome.tiff",
    "152Sm_CD66b.ome.tiff",
    "153Eu_LAG3.ome.tiff",
    "154Sm_TIM3.ome.tiff",
    "155Gd_FoxP3.ome.tiff",
    "156Gd_CD4.ome.tiff",
    "158Gd_CTLA4.ome.tiff",
    "159Tb_Brachiury.ome.tiff",
    "160Gd_CD11c.ome.tiff",
    "161Dy_CD86.ome.tiff",
    "162Dy_CD8.ome.tiff",
    "163Dy_ERG.ome.tiff",
    "164Dy_S100.ome.tiff",
    "165Ho_PD1.ome.tiff",
    "167Er_GranzymeB.ome.tiff",
    "169Tm_CD56.ome.tiff",
    "170Er_CD3.ome.tiff",
    "171Yb_CD21.ome.tiff",
    "172Yb_CD10.ome.tiff",
    "173Yb_Ki67.ome.tiff",
    "174Yb_CD72a.ome.tiff",
    "175Lu_GATA3.ome.tiff",
    "176Yb_CD80.ome.tiff",
    "191Ir_DNA1.ome.tiff",
    "193Ir_DNA2.ome.tiff",
    "208Pb_208Pb.ome.tiff",
    "209Bi_aSMA.ome.tiff",
    "80ArAr_80ArAr.ome.tiff",
]

if not roi_dir.exists():
    raise FileNotFoundError(f"ROI directory does not exist: {roi_dir}")

imgs = []
shapes = []

print(f"Reading ROI from: {roi_dir}")
print(f"Writing stacked TIFF to: {out_path}")

for fn in channel_order:
    p = roi_dir / fn
    if not p.exists():
        raise FileNotFoundError(f"Missing channel: {p}")
    img = tiff.imread(p)
    if img.ndim != 2:
        raise ValueError(f"{p} is not 2D. Shape: {img.shape}")
    imgs.append(img)
    shapes.append(img.shape)

if len(set(shapes)) != 1:
    raise ValueError(f"Channel shapes do not match: {set(shapes)}")

stack = np.stack(imgs, axis=0)  # CYX
tiff.imwrite(out_path, stack)

print(f"Wrote {out_path}")
print(f"Shape: {stack.shape}")
print(f"Dtype: {stack.dtype}")
print(f"Channels stacked: {len(channel_order)}")
```
## 6.2 Example Execution
```bash
python stack_one_roi.py ROI001_D13 data/img/ROI001_D13.tiff
```

Interpretation: The output file data/img/ROI001_D13.tiff contains all marker channels in a single structured image. This file serves as the direct input for Steinbock segmentation and downstream analysis.

---

# 7. Panel Annotation

Purpose: To define marker identity and segmentation roles for each channel in the multi-channel TIFF image.

Steinbock requires a `panel.csv` file that maps each channel to its corresponding biological marker and specifies how each channel should be used during segmentation. The order of rows in the panel file must exactly match the channel stacking order used during image construction.

## 7.1 Example panel.csv

```csv
channel,name,keep,deepcell
115In,CD20,1,
127I,127I,1,
131Xe,131Xe,1,
134Xe,134Xe,1,
138Ba,138Ba,1,
141Pr,CD68,1,
142Nd,CD138,1,
143Nd,CD47,1,2
144Nd,CD31,1,2
145Nd,Tbet,1,
146Nd,BCL2,1,
147Sm,CD44,1,2
148Nd,CD163,1,
149Sm,CD45RO,1,2
150Nd,PDL1,1,
151Eu,MHC_II,1,
152Sm,CD66b,1,
153Eu,LAG3,1,
154Sm,TIM3,1,
155Gd,FoxP3,1,
156Gd,CD4,1,
158Gd,CTLA4,1,
159Tb,Brachiury,1,
160Gd,CD11c,1,
161Dy,CD86,1,
162Dy,CD8,1,
163Dy,ERG,1,
164Dy,S100,1,
165Ho,PD1,1,
167Er,GranzymeB,1,
169Tm,CD56,1,
170Er,CD3,1,
171Yb,CD21,1,
172Yb,CD10,1,
173Yb,Ki67,1,
174Yb,CD72a,1,
175Lu,GATA3,1,
176Yb,CD80,1,
191Ir,DNA1,1,1
193Ir,DNA2,1,1
208Pb,208Pb,1,
209Bi,aSMA,1,2
80ArAr,80ArAr,1,
```
Interpretation:
- channel: Defines the metal isotope corresponding to each acquisition channel.
- name: Specifies the biological marker associated with the channel.
- keep: Indicates whether the channel should be retained for analysis (1 = keep).
- deepcell: Defines segmentation roles:
 - 1 denotes nuclear channels used for nucleus detection
 - 2 denotes membrane or cytoplasmic channels used for cell boundary delineation

Accurate panel annotation is critical, as segmentation performance depends on correct identification of nuclear and membrane markers. Incorrect assignments may result in poor segmentation quality and downstream analytical errors.

# 8. Steinbock Processing

Purpose: To perform cell segmentation, feature extraction, and data export from multi-channel IMC images.

Steinbock provides a modular command-line interface for processing Imaging Mass Cytometry data. The workflow includes segmentation, intensity quantification, morphological analysis, spatial graph construction, and export to standard data formats.

---

## 8.1 Segmentation

```bash
steinbock segment deepcell --app mesmer --minmax
```
Interpretation:
Segmentation is performed using DeepCell Mesmer, a deep learning-based model trained to identify cellular boundaries. The model integrates nuclear and membrane channels defined in the panel.csv file to generate accurate cell masks.

The output is stored in the masks/ directory as labeled TIFF images, where each pixel is assigned to a specific cell object.

## 8.2 Intensity Measurement

```bash
steinbock measure intensities
```
Interpretation:
For each segmented cell, pixel intensities are aggregated across all channels. The default aggregation function is the mean, resulting in a matrix where rows correspond to cells and columns correspond to markers.

The output is stored in the intensities/ directory.

## 8.3 Morphological Feature Extraction
```bash
steinbock measure regionprops
```
Interpretation:
Morphological properties are computed for each segmented cell, including:
- Cell area
- Centroid coordinates
- Major and minor axis lengths
- Eccentricity

These features provide structural context for downstream analyses.

The output is stored in the regionprops/ directory.

⸻

## 8.4 Spatial Neighbor Graph Construction

```bash
steinbock measure neighbors --type expansion --dmax 4
```
Interpretation:
A spatial graph is constructed where:
- Nodes represent individual cells
- Edges represent spatial proximity between cells

The parameter dmax = 4 defines the maximum expansion radius used to determine neighboring relationships.

The output is stored in the neighbors/ directory.

## 8.5 Data Export

```bash
steinbock export csv intensities regionprops -o cells.csv
```
```bash
steinbock export anndata --intensities intensities --data regionprops --neighbors neighbors -o cells.h5ad
```
```bash
steinbock export graphs --format graphml --data intensities --data regionprops
```
Interpretation:
- cells.csv: Combined table containing intensity and morphological features
- cells.h5ad: AnnData object for interoperability with Python-based single-cell workflows
- graphs/: GraphML files encoding spatial relationships between cells

These outputs form the basis for downstream analysis in R or Python environments.

# 9. Image Metadata

Purpose: To define image-level metadata required for integration of image data with cell-level measurements.

Steinbock requires an `images.csv` file that links each image to its corresponding metadata, including dimensions and identifiers. This file ensures consistency between image data, segmentation masks, and extracted features.

## 9.1 Example images.csv

```csv
image,width_px,height_px
ROI001_D13,1000,1000
```
Interpretation:
- image: Unique identifier of the ROI; must match filenames used in segmentation and intensity outputs
- width_px: Image width in pixels
- height_px: Image height in pixels

Accurate metadata alignment is critical for proper mapping between spatial coordinates, segmentation masks, and cell-level data. Mismatches in identifiers or dimensions may lead to errors during data import or spatial visualization.

## 9.2 Generation in R

```R
images <- data.frame(
  image = "ROI001_D13",
  width_px = 1000,
  height_px = 1000
)

write.csv(images, "data/images.csv", row.names = FALSE)
```
Interpretation:
The generated images.csv file enables the integration of spatial coordinates and cell-level measurements within the SpatialExperiment object used in downstream R-based analysis.

# 9. Image Metadata

Purpose: To define image-level metadata required for integration of image data with cell-level measurements.

Steinbock requires an `images.csv` file that links each image to its corresponding metadata, including dimensions and identifiers. This file ensures consistency between image data, segmentation masks, and extracted features.

## 9.1 Example images.csv

```csv
image,width_px,height_px
ROI001_D13,1000,1000
```

Interpretation:

- image: Unique identifier of the ROI; must match filenames used in segmentation and intensity outputs  
- width_px: Image width in pixels  
- height_px: Image height in pixels  

Accurate metadata alignment is critical for proper mapping between spatial coordinates, segmentation masks, and cell-level data. Mismatches in identifiers or dimensions may lead to errors during data import or spatial visualization.

## 9.2 Generation in R

```r
images <- data.frame(
  image = "ROI001_D13",
  width_px = 1000,
  height_px = 1000
)

write.csv(images, "data/images.csv", row.names = FALSE)
```

Interpretation:  
The generated `images.csv` file enables the integration of spatial coordinates and cell-level measurements within the SpatialExperiment object used in downstream R-based analysis.

---

# 10. R Analysis

Purpose: To import Steinbock-generated data into R and organize it into a structured format suitable for downstream single-cell and spatial analysis.

The `imcRtools` package provides functions to read Steinbock outputs and convert them into a `SpatialExperiment` object.

## 10.1 Load Required Libraries

```r
library(imcRtools)
library(SingleCellExperiment)
library(SpatialExperiment)
```

## 10.2 Import Steinbock Data

```r
spe <- read_steinbock(
  path = "data",
  intensities_folder = "intensities",
  regionprops_folder = "regionprops",
  graphs_folder = "neighbors",
  image_file = "images.csv",
  panel_file = "panel.csv"
)
```

Interpretation:

- assay(spe, "counts"): Raw marker intensities per cell  
- colData(spe): Cell-level metadata  
- spatialCoords(spe): Spatial coordinates  
- colPairs(spe): Cell-cell interaction graph  

## 10.3 Inspect Data Structure

```r
dim(spe)
assayNames(spe)
head(colData(spe))
rownames(spe)
```
