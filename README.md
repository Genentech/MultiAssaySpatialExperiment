# MultiAssaySpatialExperiment

Multi-Assay Experiment with spatial context for integrative analysis of spatially resolved multi-omics data.

## What is it?

**MultiAssaySpatialExperiment** extends [MultiAssayExperiment](https://bioconductor.org/packages/MultiAssayExperiment) with spatial layers: images, labels, points, and shapes. It links assay columns to spatial layer rows via an optional `spatialMap` table (`assay`, `colname`, `region`, `instance_id`), enabling multi-assay analysis in spatial context.

## Why use it?

### If you use SpatialExperiment

**SpatialExperiment** handles one assay per object. When you have *multiple assays* that share the same spatial layout (e.g., RNA + protein from the same tissue, or multiple Visium samples), you either keep separate SPE objects or manage ad-hoc lists. **MultiAssaySpatialExperiment** gives you a single container: one `ExperimentList` with multiple assays, one `colData` for sample metadata, and spatial layers (points, shapes, images) shared across assays. Subsetting by sample, assay, or spatial region propagates correctly. Coercion `as(spe, "MultiAssaySpatialExperiment")` and `as(mase, "SpatialExperiment")` keep you interoperable.

### If you use SpatialFeatureExperiment

**SpatialFeatureExperiment** extends SPE with `colGeometries`, `rowGeometries`, and `annotGeometries`—rich spatial structures per assay. When you need *multi-assay* geometry-aware analysis (e.g., gene expression and cell morphology from the same cells), MASE provides the multi-assay scaffold. Each assay can be an SFE; `spatialMap` links assay columns to shared point or shape layers. Use `annotateWithRegions` for point-in-polygon annotation and `aggregateByRegion` to aggregate assays by shape (e.g., sum expression per cell or region). Coercion to/from SFE preserves geometry when assays are SFE-compatible.

### If you use SpatialData

**SpatialData** uses a standardized on-disk schema (Zarr, Parquet) for images, labels, points, and shapes. When you want a *Bioconductor-native S4 object* for analysis (compatible with MultiAssayExperiment, SummarizedExperiment, and Bioconductor workflows), MASE fills that role. Coercion `as(sd, "MultiAssaySpatialExperiment")` and `as(mase, "SpatialData")` let you move between the SpatialData ecosystem and Bioconductor’s S4 structures. MASE adds the `ExperimentList` + `sampleMap` + `spatialMap` model that SpatialData does not provide directly.

## Key features

- **Multi-assay + spatial**: One object for multiple assays (RNA, protein, morphology, etc.) with shared points, shapes, and images
- **Instance-level mapping**: `spatialMap` links each assay column to a row in a points or shapes layer
- **Subsetting**: `subsetByColData`, `subsetByColumn`, `subsetByAssay`, `subsetByBoundingBox`, `subsetByPolygon`; spatial filters propagate to assays
- **Annotation and aggregation**: `annotateWithRegions` (point-in-polygon), `aggregateByRegion` (sum, mean, count by shape)
- **Interoperability**: Coercion to/from SpatialExperiment, SpatialFeatureExperiment, and SpatialData

## Quick example

```r
library(MultiAssaySpatialExperiment)

# Wrap a SpatialExperiment as single-assay MASE
spe <- SpatialExperiment::SpatialExperiment(...)
mase <- as(spe, "MultiAssaySpatialExperiment")

# Or build from scratch with points and shapes
mase <- MultiAssaySpatialExperiment(
    experiments = ExperimentList(rna = rna_mat, protein = protein_mat),
    colData = colData,
    sampleMap = sampleMap,
    points = PointsLayerList(centroids = coord_df),
    shapes = ShapesLayerList(cells = cell_polygons),
    spatialMap = spatialMap
)

# Point-in-polygon: which cell does each centroid belong to?
mase <- annotateWithRegions(mase, points = "centroids", shapes = "cells")

# Aggregate expression by cell
agg <- aggregateByRegion(mase, by = "cells", FUN = "sum")
```

## Documentation

See the vignettes: **Introduction**, **Subset**, and **Coercion**.
