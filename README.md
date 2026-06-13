# MultiAssaySpatialExperiment

**Multi-modal spatial transcriptomics with Bioconductor**

[![License: Apache-2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-devel-brightgreen.svg)](http://bioconductor.org/packages/devel/bioc/html/MultiAssaySpatialExperiment.html)

## Overview

**MultiAssaySpatialExperiment** extends [MultiAssayExperiment](https://bioconductor.org/packages/MultiAssayExperiment) with spatial context for integrated analysis of spatially resolved multi-omics data. It provides:

- **Multi-assay support**: Multiple experiments (RNA, protein, morphology) in one object
- **Rich spatial layers**: Images, labels (segmentation masks), points (transcripts, centroids), and shapes (cell boundaries, tissue regions)
- **Instance-level mapping**: Link assay columns to specific spatial features via `spatialMap`
- **Built-in readers**: Load data from Xenium, Visium, Visium HD, MERSCOPE, and CosMx
- **Spatial operations**: Annotate, aggregate, subset by spatial criteria
- **Full interoperability**: Coercion to/from SpatialExperiment and SpatialFeatureExperiment

## Installation

```r
# From Bioconductor (devel)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MultiAssaySpatialExperiment", version = "devel")

# Or install from a local source checkout
# R CMD INSTALL path/to/MultiAssaySpatialExperiment
```

## Quick Start

### Load vendor data directly

```r
library(MultiAssaySpatialExperiment)

# Read 10x Xenium data
mase <- readXeniumMASE("path/to/xenium/output")

# Read 10x Visium HD data
mase <- readVisiumHDMASE("path/to/visium_hd/output",
                         bin_size = c("008", "016"))

# Read Vizgen MERSCOPE data
mase <- readMERSCOPEMASE("path/to/merscope/output",
                         segmentation = "cellpose",
                         load_transcripts = TRUE)

# Read NanoString CosMx data
mase <- readCosMxMASE("path/to/cosmx/output")

# Read standard 10x Visium data
mase <- readVisiumMASE("path/to/visium/output")
```

**Supported technologies**: Xenium, Visium, Visium HD, MERSCOPE (Vizgen), CosMx (NanoString)

### Build from components

```r
# Create from scratch
mase <- MultiAssaySpatialExperiment(
    experiments = ExperimentList(
        rna = rna_counts,
        protein = protein_counts
    ),
    colData = sample_metadata,
    sampleMap = sample_map,
    points = PointsLayerList(
        transcripts = transcript_coords,
        centroids = cell_centroids
    ),
    shapes = ShapesLayerList(
        cells = cell_polygons,
        tissue = tissue_boundary
    ),
    spatialMap = spatial_map
)
```

### Convert from existing objects

```r
# From SpatialExperiment
spe <- SpatialExperiment(...)
mase <- as(spe, "MultiAssaySpatialExperiment")

# From SpatialFeatureExperiment
sfe <- SpatialFeatureExperiment(...)
mase <- as(sfe, "MultiAssaySpatialExperiment")
```

## Core Features

### Spatial Operations

```r
# Point-in-polygon annotation
mase <- annotateWithRegions(mase,
                            points = "centroids",
                            shapes = "cells")

# Aggregate expression by spatial region
cell_expr <- aggregateByRegion(mase,
                               by = "cells",
                               FUN = "sum")

# Spatial subsetting
library(sf)
roi <- st_polygon(list(matrix(c(0, 0, 1000, 0, 1000, 1000, 0, 1000, 0, 0),
                              ncol = 2, byrow = TRUE)))
mase_tissue <- subsetByPolygon(mase, roi)

mase_bbox <- subsetByBoundingBox(mase,
                                 xmin = 0, xmax = 1000,
                                 ymin = 0, ymax = 1000)
```

### Standard Subsetting

Uses the same `subsetBy*` API as MultiAssayExperiment — no separate spatial query API.
Specimen subsetting (`[, j]`, `subsetByColData`) and observation filtering
(`subsetByColumn` with a list) both propagate to `spatialMap`, `imgData`, and linked
points, shapes, images, and labels.

```r
# By sample, assay, column, or observation metadata
mase[, colData(mase)$tissue_type == "tumor"]
mase[c("rna", "protein"), ]
subsetByAssay(mase, "rna")
cdf <- colData(experiments(mase)[["rna"]])
subsetByColumn(mase, list(rna = cdf$region == "core"))
```

### Construction helpers

```r
spmap <- buildSpatialMap(sampleMap(mase), region = "cells", element_type = "shapes")
prepared <- prepMASE(experiments(mase), colData(mase), sampleMap(mase),
                     points = spatialPoints(mase), spatialMap = spmap)
mase2 <- do.call(MultiAssaySpatialExperiment, prepared)
```

### Access Spatial Data

```r
# Points (DataFrame with coordinates)
transcripts <- spatialPoints(mase)$transcripts
centroids <- spatialPoints(mase)$centroids

# Shapes (sf objects with geometries)
cell_boundaries <- spatialShapes(mase)$cells
tissue_outline <- spatialShapes(mase)$tissue

# Images and labels
spatialImages(mase)
spatialLabels(mase)

# Spatial mapping table
spatialMap(mase)
```

## Why MultiAssaySpatialExperiment?

### For SpatialExperiment users

**Problem**: SpatialExperiment handles one assay per object. Multi-assay data requires managing separate objects or ad-hoc lists.

**Solution**: MASE provides a unified multi-assay container with shared spatial layers. Subset once, filter everywhere. Full coercion support maintains interoperability.

```r
# Instead of managing multiple SPE objects
spe_rna <- SpatialExperiment(...)
spe_protein <- SpatialExperiment(...)

# Use one MASE object
mase <- MultiAssaySpatialExperiment(
    experiments = ExperimentList(rna = spe_rna, protein = spe_protein),
    ...
)
```

### For SpatialFeatureExperiment users

**Problem**: SFE extends SPE with rich geometry support but is single-assay.

**Solution**: MASE provides the multi-assay scaffold. Each assay can be an SFE. Instance-level `spatialMap` links columns across assays to shared spatial features.

```r
# Multi-assay with SFE-compatible structure
mase <- MultiAssaySpatialExperiment(
    experiments = ExperimentList(
        rna = sfe_rna,        # SFE object
        morphology = sfe_morph # SFE object
    ),
    shapes = ShapesLayerList(cells = boundaries),
    spatialMap = map  # Links both assays to same cells
)
```

## Architecture

### Core Components

**Inherited from MultiAssayExperiment**:
- `ExperimentList`: Multiple experiments (SummarizedExperiment, etc.)
- `colData`: Sample-level metadata
- `sampleMap`: Links experiment columns to samples

**Spatial extensions**:
- `points`: PointsLayerList (transcripts, centroids, etc.)
- `shapes`: ShapesLayerList (cells, nuclei, regions)
- `images`: RasterLayerList (user-attached rasters; distinct from reader metadata)
- `labels`: RasterLayerList (segmentation masks)
- `imgData`: Specimen-level image metadata from vendor readers (SPE-compatible)
- `spatialMap`: Instance-level mapping (assay/colname → spatial feature)

### The spatialMap Table

Links assay columns to spatial features at instance level:

| assay | colname | element_type | region | instance_id |
|-------|---------|--------------|--------|-------------|
| rna   | ACGT-1  | shapes       | cells  | cell_001    |
| rna   | ACGT-2  | shapes       | cells  | cell_002    |
| protein | ACGT-1 | shapes      | cells  | cell_001    |

- **Foreign key** from (assay, colname) to (element_type, region, instance_id)
- Enables multi-assay analysis on same spatial features
- Optional but powerful for integrated analysis

## Data Readers

Built-in support for major spatial transcriptomics platforms:

| Technology | Reader Function | Data Types | Special Features |
|------------|----------------|------------|------------------|
| **10x Xenium** | `readXeniumMASE()` | Counts, cells, boundaries, transcripts | Auto-format detection (HDF5/Parquet) |
| **10x Visium** | `readVisiumMASE()` | Counts, positions, images | H&E integration, spot geometries |
| **10x Visium HD** | `readVisiumHDMASE()` | Multi-bin counts, images | Multiple resolutions (002µm, 008µm, 016µm) |
| **Vizgen MERSCOPE** | `readMERSCOPEMASE()` | FOV-based, transcripts, boundaries | Multi-FOV, cellpose/watershed segmentation |
| **NanoString CosMx** | `readCosMxMASE()` | FOV-based, expression, boundaries | Multi-FOV, GeoParquet geometries |

All readers:
- Return fully-validated MASE objects
- Handle multiple file formats (HDF5, Parquet, CSV, GeoJSON)
- Support optional components (transcripts, images, labels)
- Use S4 generics (extensible for database backends)

## Documentation

### Vignettes

1. **Introduction to MultiAssaySpatialExperiment**: Overview, construction, basic operations
2. **Working with MultiAssaySpatialExperiment**: Subsetting, spatial annotation, aggregation, labels ↔ shapes
3. **MultiAssaySpatialExperiment Use Cases**: Real data workflows with readers, multi-assay integration, coordinate transforms
4. **MultiAssaySpatialExperiment Cheatsheet**: Quick reference for common tasks

```r
browseVignettes("MultiAssaySpatialExperiment")
```

### Key Functions

**Data import**:
- `readXeniumMASE()`, `readVisiumMASE()`, `readVisiumHDMASE()`
- `readMERSCOPEMASE()`, `readCosMxMASE()`

**Spatial operations**:
- `annotateWithRegions()`: Point-in-polygon annotation
- `aggregateByRegion()`: Aggregate assays by spatial features
- `subsetByBoundingBox()`: Subset to spatial extent
- `subsetByPolygon()`: Subset to polygon interior
- `subsetByColumn()`: Filter assay columns (with spatial propagation when `y` is a list)
- `spatialJoin()`: Join spatial layer tables (DataFrame × DataFrame via `sf`)

**Accessors**:
- `spatialPoints()`, `spatialShapes()`, `spatialImages()`, `spatialLabels()`
- `spatialMap()`, `imgData()`

**Coercion**:
- `as(x, "MultiAssaySpatialExperiment")`: From SPE/SFE
- `as(mase, "SpatialExperiment")`: To SPE

## Design Principles

### Terminology and spatialdata alignment

MASE aligns with [spatialdata](https://spatialdata.scverse.org/) (Python) for interoperability:

- **Specimen**: a row in `colData` (primary identifier); often a tissue section,
  patient, or replicate
- **Observation**: a column in an assay matrix (cell, spot, bin)
- **LayerList** vs **Element**: MASE uses "LayerList" (PointsLayerList, ShapesLayerList) to avoid R name collisions. Maps to spatialdata "elements" conceptually.
- **element_type**: Discriminates point vs shape layers in `spatialMap` (spatialdata uses slot names; MASE needs explicit column)
- **region**: Layer name within element type (spatialdata "region", MASE "layer" — same concept)
- **instance_id**: Row identifier within layer (consistent across both)

See `?MultiAssaySpatialExperiment` for detailed terminology documentation.

### Instance-Level Validation

`instance_id` validation ensures referential integrity:
- Must exist in target spatial layer
- Cannot contain NAs
- Recommended types: character, integer, factor
- Informational warnings for unusual types (numeric, complex)

Helps catch data inconsistencies early.

## Performance and Extensibility

### S4 Generics for File I/O

All readers use S4 generics (`readParquetForMASE()`, `readGeoParquetForMASE()`, etc.) enabling:

- **Package extension**: Other packages can register optimized methods
- **Consistent interface**: Technology readers use same primitives

### Component-Based Architecture

Technology readers orchestrate reusable components:
- **Component readers**: Format-agnostic (Parquet, CSV, HDF5, GeoJSON)
- **Technology orchestrators**: Handle platform-specific logic
- **Configuration registry**: Extensible for new technologies

Add new platforms without modifying core code.

## Package Statistics

- **20 R source files**
- **11 test files**
- **4 vignettes**
- **5 technology readers**
- **Coercion**: SpatialExperiment and SpatialFeatureExperiment ↔ MASE

## Development Status

**Version 0.9.3** - Pre-release for community feedback

- ✅ Core class implementation complete
- ✅ Spatial operations (annotate, aggregate, subset)
- ✅ Reader architecture (5 technologies)
- ✅ Full interoperability (SPE, SFE)
- ✅ Comprehensive documentation
- ⏳ Bioconductor submission planned

## Getting Help

- **Documentation**: `?MultiAssaySpatialExperiment`
- **Vignettes**: `browseVignettes("MultiAssaySpatialExperiment")`
- **Issues**: [GitLab Issues](https://code.roche.com/GP/MultiAssaySpatialExperiment/-/issues)

## License

Apache License 2.0 — see DESCRIPTION.

Copyright (c) 2023-2026 Genentech, Inc.

## Related Projects

- [MultiAssayExperiment](https://bioconductor.org/packages/MultiAssayExperiment): Parent class for multi-assay data
- [SpatialExperiment](https://bioconductor.org/packages/SpatialExperiment): Single-assay spatial data
- [SpatialFeatureExperiment](https://bioconductor.org/packages/SpatialFeatureExperiment): Geometry-aware single-assay
