# MultiAssaySpatialExperiment 0.9.6

## Documentation

- Reconciled the vignettes into a numbered, user-first set:
  - *1. Introduction to MultiAssaySpatialExperiment* --- overview, anatomy,
    construction, accessors, subsetting, spatial operations, coercion, and a
    compact quick-reference table.
  - *2. Working with MultiAssaySpatialExperiment* --- construction patterns,
    subsetting, spatial annotation and aggregation, and labels <-> shapes.
  - *3. MultiAssaySpatialExperiment use cases* --- platform readers (Xenium,
    Visium, Visium HD, CosMx, MERSCOPE) and real-world spatial workflows.
- Removed the *Cheatsheet* vignette (its chunks were globally unevaluated, which
  `BiocCheck` flags); its quick-reference content is now a table in the
  Introduction.
- Applied Bioconductor vignette conventions: numbered `VignetteIndexEntry`,
  `BiocStyle::doc_date()` compiled dates, and structured author metadata.
