# MultiAssaySpatialExperiment 0.99.8

## Internal changes

- `subsetByBoundingBox()` and `subsetByPolygon()` were structurally identical
  apart from the functions used to index each layer, and each contained four
  near-duplicate loops. Both now delegate to one `.subsetMASEByGeometry()`
  engine, which also removes redundant work: each method previously indexed
  every spatial layer **twice**, once to build the retained-instance map and
  again, with identical arguments, to build the filtered layers, so
  `subsetByPolygon()` ran its point-in-polygon join twice per layer on every
  call. Equivalence was checked against the previous implementation across
  partial-hit, everything-retained and empty-result cases for both methods,
  comparing dimensions, layer contents, layer names and `spatialMap`.
- The review asked whether the package's `for` loops could use the apply
  family. Eight are removed by the refactor above. Of the remainder, fifteen
  are file-pattern searches that `break` on the first hit, where an apply
  variant cannot short-circuit and would test every candidate path; the rest
  either fold state across iterations or exist for their side effects, where a
  loop is the clearer form. The two examples named in the review are kept:
  `.get_img_idx_mase()` is a near-verbatim port of
  `SpatialExperiment:::.get_img_idx()`, deliberately structured like its
  upstream counterpart, and `.import_SFE_generics()` is a two-line side-effect
  loop registering generics.

# MultiAssaySpatialExperiment 0.99.7

## Documentation

- Expanded the narrative in *MultiAssaySpatialExperiment use cases*, which the
  Bioconductor review found short on explanation of what each step does and
  why. The reader sections now describe the arguments that actually change the
  result and what they cost; the coercion sections say why you would coerce
  rather than re-read the vendor output, and what survives each conversion;
  Workflow 3 states the problem it solves instead of listing what could be
  done; Workflow 4 explains why its steps are in that order, since the point of
  that section is the ordering; and the summary says plainly what the container
  does not do, which is also why the affine-transformation example is
  deliberately plain R.
- Removed three further blocks of explanatory comments from code chunks
  (`# The conversion:` twice and `# The combined MASE has:`), matching the ones
  already removed in 0.99.7. No `cat()` call or comment block now stands in for
  prose in any vignette.
- Copy edit across all five vignettes for correctness, consistency, clarity and
  concision. Substantive fixes: three cross-references named a vignette
  *MultiAssaySpatialExperiment Use Cases* when the title is lower case; one
  sentence in the introduction said an assay column belongs to a "sample" two
  lines after establishing that it belongs to a *specimen*, which is the exact
  distinction that section teaches. Consistency fixes: vignette
  cross-references were split between bold and italic and are now uniformly
  italic, as elsewhere in the suite; the two "See also" lists disagreed on
  heading, ordering and capitalisation and are now identical in form; two
  Oxford commas were brought into line with the prevailing style. One 60-word
  sentence was split.
- Restructured the vignettes in response to the Bioconductor review, which
  observed that the introduction read like a developer vignette. The
  introduction is now motivation followed by a single worked example carried
  from construction through subsetting, annotation and aggregation to
  coercion, matching the shape of `SingleCellExperiment`'s `intro.Rmd`. The
  slot-by-slot anatomy, the mapping tables and the entity-relationship
  diagram move to a new *Design of MultiAssaySpatialExperiment* vignette,
  which is both the split `SingleCellExperiment` makes between `intro` and
  `devel` and the `<Package>-design` convention already used by the sibling
  packages in this suite.
- Added a *MultiAssaySpatialExperiment cheatsheet* vignette, a one-page API
  reference modelled on `MultiAssayExperiment`'s. Two vignettes previously
  pointed readers at a cheatsheet that did not exist.
- Vignette files follow the suite naming convention, so pkgdown links the
  introduction as the package landing page:
  `Introduction.Rmd` becomes `MultiAssaySpatialExperiment.Rmd`,
  `WorkingWithMASE.Rmd` becomes `MultiAssaySpatialExperiment-workflows.Rmd`,
  and `UseCases.Rmd` becomes `MultiAssaySpatialExperiment-usecases.Rmd`.
- The entity-relationship diagram is rebuilt at roughly four times the
  resolution, with entities declared in the order the text introduces them and
  genuine per-entity colouring. The prose had described "tan" and "blue" boxes,
  but the source set a single colour and defined no per-entity styling, so no
  such distinction was ever rendered. `spatialImages` and `spatialLabels` are
  included for the first time; their absence was part of why the diagram and
  its description could not be reconciled. Crow's-foot notation is now
  explained.
- Added a prerequisites section stating what background is assumed, and that
  `sf` is required only for the `shapes` slot.
- Replaced the remaining bullet lists that stood in for explanation with
  prose, removed the `cat()` calls that printed static workflow text from code
  chunks, gave the subsetting, spatial-operation and coercion sections
  runnable examples, and renamed headings that merely repeated the sentence
  below them.
- The subsetting examples use two assays and distinguish specimen identifiers
  (`P1`, `P2`) from observation identifiers (`S1`...`S4`). With one assay whose
  column names matched its specimen names, subsetting by specimen and by assay
  column produced identical output and the distinction was invisible.
- Fixed bullet lists that were not separated from the preceding line, which
  pkgdown did not render as lists.
- Renamed "Next steps" sections to "See also", as they enumerate the other
  vignettes rather than a progression.

## Bug fixes

- Added `Spatial` to `biocViews`.
- `.check_visium_files()` no longer builds a message with `paste0()` inside
  the condition signal; the subdirectory name is computed once and reused.

## Internal changes

- Removed four unused internal functions: `.spatialLabels_from_rasters()`,
  `.read_label_rasters()`, `set_block_size()` and `list_technologies()`.
- Fixed the indentation of `else` continuations in `.mergeMASE()`, the `cbind`
  method, and the `SpatialFeatureExperiment` coercion.

# MultiAssaySpatialExperiment 0.99.6

## Bug fixes

- Fixed a spurious `no definition for class "MultiAssaySpatialExperiment"`
  warning (one per `setAs()`/`setMethod()` call registered in
  `.register_SFE_coercion()`/`.register_SFE_methods()`, ~20 total) emitted
  whenever `SpatialFeatureExperiment` was already loaded before
  `MultiAssaySpatialExperiment` (e.g. `library(SpatialFeatureExperiment);
  library(MultiAssaySpatialExperiment)`). Registration in that case runs
  synchronously inside `.onLoad`, which always completes *before* this
  package is attached to `search()`; `.sfe_ext_env` (added in 0.99.5 to
  work around namespace sealing) was parented to `globalenv()`, so at
  registration time it couldn't yet see `"package:MultiAssaySpatialExperiment"`
  to validate the class in the method signature. `.sfe_ext_env` is now
  parented to this package's own namespace instead, where the class is
  always immediately reachable regardless of attachment timing. Harmless
  (registration and dispatch both still worked correctly), but noisy.

# MultiAssaySpatialExperiment 0.99.5

## Bug fixes

- `as(spe/sfe, "MultiAssaySpatialExperiment")` now store the assay as a plain
  `SingleCellExperiment` and keep all spatial content (colGeometries,
  annotGeometries, `unit`, col-margin `spatialGraphs`) in MASE's own
  `points`/`shapes`/`spatialMap`/`colPair()` layers, instead of retaining the
  full `SpatialExperiment`/`SpatialFeatureExperiment` object. Previously,
  `writeParquet()` (BiocDuckDB) only ever read those MASE-native layers, so
  anything living only in a nested SFE's colGeometries beyond the primary
  centroids, annotGeometries, or spatialGraphs was silently dropped on
  write even though it looked preserved in memory. `as(mase,
  "SpatialExperiment"/"SpatialFeatureExperiment")` and the SFE-delegating
  methods (`dimGeometry`, `annotGeometry`, `bbox`, `unit`, `splitByCol`,
  `aggregate`, `findDebrisCells`) now reconstruct a real SPE/SFE view from
  those native layers on demand; `spatialGraphs`/`spatialGraph`/
  `spatialGraphNames` read `colPair()` directly. `localResults` is not
  carried over by this coercion (no MASE-native home for it yet); it now
  errors with a message that says so instead of a misleading
  "no SpatialFeatureExperiment assay" message.
- Fixed a registration bug where the entire `SpatialFeatureExperiment`
  integration (all coercions and the `dimGeometry`/`unit`/`bbox`/etc.
  delegating methods) would silently fail to register whenever anything
  loaded `MultiAssaySpatialExperiment` before the user loaded
  `SpatialFeatureExperiment` in the same session (for example, loading
  `BiocDuckDB` first, since it imports this package) -- `as(sfe,
  "MultiAssaySpatialExperiment")` would then silently fall back to the
  generic `SpatialExperiment` coercion with no error. The deferred
  registration hook ran after this namespace was already sealed, and
  `setAs()`/`setMethod()` need to create a new top-level binding the first
  time a generic is used from a namespace; that's forbidden once a
  namespace is sealed, and the resulting error was masked by an
  overly-broad `grepl("locked", ...)` catch in the load hook. Registration
  now targets a small dedicated, never-sealed environment instead of this
  namespace, which works regardless of load order; the test suite's own
  `pkgload::load_all()` workaround for this (`helper-sfe-register.R`) is no
  longer needed and has been removed.

## Documentation

- Added the missing `@aliases coerce,...-method` entries for the
  `list`/`List` -> `RasterLayerList`/`PointsLayerList`/`ShapesLayerList`
  coercions in `?SpatialLayerList-class`; these were implemented but
  undocumented, which `R CMD check` flagged as missing documentation
  entries.

# MultiAssaySpatialExperiment 0.99.4

## Documentation

- Replaced em dashes with commas or colons in the vignettes.

# MultiAssaySpatialExperiment 0.99.3

## Documentation

- Added `URL` and `BugReports` fields to DESCRIPTION.
- Added a package-level man page (`?MultiAssaySpatialExperiment`).

# MultiAssaySpatialExperiment 0.99.2

## Documentation

- Added `\value` sections (roxygen `@return`) to the `MultiAssaySpatialExperiment`
  class, combine, subset, and `SpatialExperiment`-coercion man pages, and to
  `SpatialLayerList-class`, documenting the values returned by the constructors,
  accessors, replacement, coercion, combining, and subsetting methods. Resolves
  the `R CMD BiocCheck` "missing \value" WARNING.

# MultiAssaySpatialExperiment 0.99.1

## Bug fixes

- `readMERSCOPEMASE()` and `readCosMxMASE()` now have runnable examples backed
  by minimal mock output directories bundled under `inst/extdata/`
  (`merscope_mock/`, `cosmx_mock/`; generated by `inst/scripts/make-example-mocks.R`).
  This raises the share of exported-object man pages with runnable examples above
  the 80% BiocCheck threshold, resolving the `R CMD BiocCheck` ERROR.

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
