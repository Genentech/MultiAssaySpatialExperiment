# MultiAssaySpatialExperiment 0.99.16

## Internal changes

- Converting a `SpatialExperiment` or a `SpatialFeatureExperiment` adds a helper
  column to a copy of the sample map, which `buildSpatialMap()` reads the instance
  ids from. It is named `__instance_id__`. It used to be `.instance_id`, and now
  follows the style of the columns that reach files and schemas, such as
  `__name__` and `__index__`. Only the name changed: the helper column is not
  kept, and the `instance_id` column of the spatial map is the same.

# MultiAssaySpatialExperiment 0.99.15

## Bug fixes

- Points converted from a `SpatialData` had no `x` or `y` column: the
  coordinates stayed in spatialdataR's point geometry, next to the Python index
  column, so the layer could not be used by MASE's own operations. For
  example, `subsetByBoundingBox()` over the whole extent of the blobs example
  took its points layer from 200 rows to none. The geometry is now turned into
  `x` and `y` (and `z` when there is one), placed first, and every other
  column, including `__null_dask_index__` and any feature column, is kept as
  stored, so the element still goes back unchanged.
- `c()` of two `MultiAssaySpatialExperiment` objects lost the coordinate
  transforms of the second one. Joining the two metadata lists left two
  entries called `transforms`, and only the first is ever found. The
  transforms of both objects are now kept in one entry. They are filed per
  element and element names are unique across the two objects, so they do not
  collide.
- Renaming an element with `pointNames<-`, `shapeNames<-`, `imageNames<-` or
  `labelNames<-` left the element's coordinate transforms filed under its old
  name, where nothing found them (for points and shapes it already updated the
  `spatialMap`). The transforms now follow the rename.

# MultiAssaySpatialExperiment 0.99.14

## Bug fixes

- A `SpatialData` points or shapes element now links to its table through the
  table's `instance_key`, whatever that column is called. Before, only an
  `instance_id` or `id` column worked, so a table keyed on, say, `cell_id`
  produced a `spatialMap` that pointed at nothing.
- Converting to `SpatialData` could give a shapes element with no rows, with
  no error: a table annotating the shapes was read and written with its region
  and instance ids in `colData`, where spatialdataR keeps them in `int_colData`
  and matches them to the element through the element's `instance_key`. spatialdataR then
  kept none of the shapes. Tables now carry their keys where spatialdataR looks,
  elements declare their `instance_key`, and a table built by spatialdataR is
  read the same way. A layer that loses rows in the `SpatialData` now warns.
- The CRS of points and shapes was dropped, silently, going to `SpatialData`. It
  is now carried over.
- A failed conversion to `SpatialData` removed none of the DuckDB tables it had
  written. It now does. (After a successful one the tables stay on the
  connection, because the `SpatialData` reads from them.)
- Unloading and loading `SpatialFeatureExperiment` or `spatialdataR` again left
  the integration's methods missing without a word, because the record that
  they were registered outlived the old namespace. Registration is now tied to
  the namespace that is loaded.
- The `spatialGraphs()`, `spatialGraph()` and `spatialGraphNames()` methods
  ignored `sample_id` and returned one merged graph. They now read one sample at
  a time and return the shape SpatialFeatureExperiment does.
- The spatialdataR accessors agree with spatialdataR: looking up a name or
  position that is not there is an error, `NULL` empties a layer, and
  `element(x, i) <- value` takes a position. `point(x, 10) <- v` no longer
  invents an element called `point3`. The names of images, labels, points and
  shapes can be replaced, and the `spatialMap` follows.
- Going to `SpatialData` now warns that the sample-level `colData` is dropped,
  and a layer whose instance ids do not match the assay's column names warns
  that its rows are taken in stored order.
- `dimGeometries()` no longer declares a default for `MARGIN` that could never
  apply.
- The metadata of a `SpatialFeatureExperiment` was dropped by the coercion, and
  so was the weight style of every sample's graph but the first. Both now
  survive the round trip.
- `findDebrisCells()` on a `MultiAssaySpatialExperiment` now returns a
  `MultiAssaySpatialExperiment` with the `is_debris` flag in the experiment's
  `colData`, as `aggregate()` and `splitByCol()` already did. It returned a
  bare `SpatialFeatureExperiment` before.
- `localResults()` agrees with `localResultNames()`: with nothing requested it
  returns an empty list instead of an error. Asking for a named result still
  fails, with a message that says local results are not carried over.

- Reloading the package no longer stacks duplicate load hooks for
  SpatialFeatureExperiment and spatialdataR, and unloading removes them
  (new `.onUnload`). A registration that fails now warns, naming the
  integration and the cause, instead of printing a bare hook error.
- The SpatialData integration never activated against the real package. It
  waited for a package called `SpatialData`, but the Bioconductor package is
  `spatialdataR`, and it called constructors (`PointFrame()`, `ShapeFrame()`)
  that spatialdataR no longer has. It now targets `spatialdataR`:
  - Points and shapes are materialized with `as.data.frame()`, because
    spatialdataR holds them as lazy DuckDB queries, and are built with
    `SpatialDataPoint()` and `SpatialDataShape()`.
  - `labels()` is registered as `base::labels(object, ...)`, which is what it
    is in spatialdataR.
  - `layer(x, i)` returns the layer that holds element `i`, and `element(x, i)`
    takes an element name or a position counted across the layers, as in
    spatialdataR. The earlier `layer(x, "images")` and
    `element(x, "images", name)` forms are gone. `element<-` is new.
  - The registration now works whichever package loads first. Before, it only
    worked if spatialdataR was loaded before this package, because it tried to
    define objects in an already sealed namespace.
  - Images and labels are still not converted. Both directions now say so in a
    warning that names them, and a table that annotates an element that was
    not converted has that region annotation removed when converting back,
    instead of producing a SpatialData that spatialdataR rejects.
  - The integration used `int_metadata()` without importing it, and called
    `layer()` and `element()` without qualifying them. It only worked in a
    session that had attached SingleCellExperiment and spatialdataR. Elsewhere
    a table that annotates a points or shapes element came out with no
    `spatialMap` at all, without an error, because one call sat inside a
    `tryCatch()` that returned `NULL`, and `element()` failed outright when
    spatialdataR was loaded but not attached.
  - Converting to SpatialData with more than one points or shapes element gave
    every element the data of the last one built, because spatialdataR's
    constructors write all elements of a kind to one table. Each element is now
    written to a table of its own, and its rows are checked afterwards. Element
    names must now be unique across points, shapes and tables, and a clash is
    reported instead of failing with an unrelated message.
  - Points whose coordinate columns are not `x` and `y`, such as the `imageX`
    and `imageY` of a MASE made from an SFE, now convert; they use the first two
    columns other than `instance_id`, as the rest of the package does.
  - A `spatialMap` that covers only some of an assay's columns no longer throws
    away the instance ids of the rest: the mapped columns keep theirs and the
    others get `NA`, with a warning. A table that annotates several regions is
    read and written per column. A table keeps its own region and instance keys,
    and column names such as `__null_dask_index__` are no longer made
    syntactic.
  - `images<-`, `labels<-`, `points<-`, `shapes<-` and `tables<-` accept the
    layer lists their getters return, not only a plain `list`.
  - Coordinate transformations of points and shapes are not converted; a
    warning now names the elements whose transformations are not the identity.
  - The integration's code moved from `inst/scripts`, where it was sourced at
    run time and so never seen by `R CMD check` or BiocCheck, into
    `R/MultiAssaySpatialExperiment-SpatialData.R`, laid out like the
    SpatialFeatureExperiment one: its generics are imported into a private
    environment and every method is registered there with an explicit `where`.
    It no longer shadows `setMethod()`, `setReplaceMethod()` and `setAs()`.
  - `spatialdataR` is now in `Suggests`, and the integration's tests moved from
    `inst/scripts` to `tests/testthat`, where `R CMD check` runs them. They
    include two clean R sessions, one per load order, which is the only way to
    see the failures above.

- Converting a `SpatialFeatureExperiment` to a MASE had several problems that
  showed up on anything but a single sample with simple graphs:
  - With several samples, `colGraph()` and `colGeometry()` were called with
    their defaults, which return only the first sample, and the coercion
    failed or was silently partial. Objects with several samples now convert
    both ways: geometries are taken for all samples, and each sample's column
    graph is merged into one graph over all columns and split again per sample.
  - A graph with the `"raw"` weight style made the converted MASE unusable
    (`"Style raw invalid"`), because the graph was rebuilt with `nb2listw()`,
    which does not accept it. It is now assembled directly, as SFE does.
  - `unit()` rebuilt a whole `SpatialFeatureExperiment` on every call and
    returned `NULL` if that failed. It now reads the unit from the object's
    metadata, or from the assay when that is itself an SFE.
  - Row geometries, row graphs and annotation graphs were dropped without a word.
    A warning now names them, as the message about `localResults` already did.

## Documentation

- The help pages for both integrations list their methods as aliases, so
  `?spatialGraphs` and the like find them, and say where the results differ from
  the upstream objects (`dimGeometry()` columns, `bbox()`, graphs, DuckDB tables).
  The README no longer carries file counts and points to the right issue
  tracker, and the use-cases vignette no longer loads `ggplot2`, which it did
  not use.
- New help page `?MultiAssaySpatialExperiment-SpatialFeatureExperiment`, which
  error messages already pointed to. It covers what the coercion keeps and
  drops, what each method does, the cost of rebuilding a view per call, and how
  `unit()` and `localResults()` differ from SFE's.
- The docs say the `spatialdataR` and `SpatialFeatureExperiment` integrations
  need the package *loaded*, not just installed, since they register when it
  loads. The package description now mentions `spatialdataR`, and the README
  no longer carries a hard-coded version.

- Documented the spatialdataR integration: a new help page,
  `?"coerce-SpatialData"`, a "From SpatialData" section in the use-cases
  vignette, two rows in the cheatsheet's coercion table, and a README section.
  The README's coercion list also now includes `as(mase,
  "SpatialFeatureExperiment")`, which was missing.
- Corrected the terminology sections that said MASE has no bare accessors, and
  that called a named entry in a list a "layer". As in spatialdataR, a layer is
  one of images, labels, points, shapes or tables, and an element is a named
  entry in a layer.

# MultiAssaySpatialExperiment 0.99.13

## Bug fixes

- `as(mase, "SpatialFeatureExperiment")` failed with `'arg' should be one of
  "full_res_image_pixel", "micron"` for any MASE built from an SFE that
  stores its unit as `"full_res_image_pixels"`, which is how older SFE
  objects, including the `SFEData` example datasets, record it. The SFE
  constructor only accepts the singular spelling, and the coercion passed the
  stored value through unchanged. The legacy spelling is now normalized both
  when the unit is stored in `metadata(mase)` and when it is read back, so
  MASE objects written before this fix also convert.

# MultiAssaySpatialExperiment 0.99.12

## Documentation

- Replaced `\dontrun{}` with `\donttest{}` in the four platform-reader
  examples that need real vendor output (`readHDF5ForMASE()`,
  `readVisiumMASE()`, `readVisiumHDMASE()`, `readXeniumMASE()`), per
  `BiocCheck`: `\dontrun{}` examples never run, even under
  `R CMD check --run-donttest` or interactive `example()`, while
  `\donttest{}` ones do, which is what lets Bioconductor's periodic checks
  catch these examples rotting out from under a real file-format change.
  Verified directly with `BiocCheck`: no ERROR or WARNING, and its
  dontrun/donttest usage NOTE now reports only `\donttest{}` tags.

# MultiAssaySpatialExperiment 0.99.11

## Documentation

- A fresh review of the current vignettes against the standards the
  Bioconductor review has consistently applied (prose over bullets, runnable
  code over listings, consistent heading conventions) found four residual
  issues, all fixed:
  - The introduction vignette's own "Next steps" section was never actually
    renamed to "See also", unlike its two siblings: it still read
    `# Where to go next`, in a full-sentence bullet style that didn't match
    the "*Title*: short phrase" convention `-workflows.Rmd` and
    `-usecases.Rmd` already use consistently. Brought into line.
  - *Working with MultiAssaySpatialExperiment*'s "Validation checklist"
    section listed construction invariants but never demonstrated any of
    them, the same "listing instead of demonstrating" pattern already fixed
    elsewhere in that vignette. Added two runnable, verified examples: a
    `sampleMap` mismatch that is silently dropped with a message during
    harmonization, and a broken `spatialMap` that fails `validObject()`
    outright, since a dangling geometry reference cannot be silently
    dropped the way a sampleMap row can.
  - Converted six narrative/comparison bulleted lists in `-usecases.Rmd` and
    `-workflows.Rmd` to prose (transcript use cases, two Xenium/Visium
    background comparisons, a data-preparation listing, the aggregation
    strategy rationale, and the rasterization/vectorization preview).
    Structural lists (vendor output file listings, parameter tables, the
    numbered platform-integration rationale, "See also") were left alone;
    those are legitimate uses of a list, not the pattern being fixed.
  - The entity-relationship diagram's `experiments(mase)` box had no visible
    connection to the `ExperimentList` class named in the surrounding prose;
    relabeled to `experiments(mase) (ExperimentList)`, matching the
    `Experiment (SummarizedExperiment)` convention the same diagram already
    uses for a different entity.
- While regenerating the diagram above, found that `build_diagrams.sh`
  itself would silently regenerate it at 784x527, the exact low resolution
  from a previous round of review, since the script never recorded the
  scale factor actually used to produce the checked-in 3168x2128 image, and
  its default background (`transparent`) rendered incorrectly in at least
  one viewer. Both fixed (`-s 4`, default background changed to `white`);
  the diagram is now reproducible from the script as checked in.

# MultiAssaySpatialExperiment 0.99.10

## Documentation

- Fixed the `R CMD check` "Rd cross-references" NOTE: 34 roxygen `\linkS4class{}`
  references to classes defined in other packages (`DataFrame` and `List`/
  `SimpleList` in `S4Vectors`; `SpatialExperiment`; `SingleCellExperiment`;
  `MultiAssayExperiment` and `ExperimentList` in `MultiAssayExperiment`;
  `SummarizedExperiment`) had no package anchor, across 7 source files and
  the 12 man pages they generate. `\linkS4class{X}` only resolves without a
  warning for classes in the current package; cross-package references now
  use `\link[pkg:X-class]{X}`. Classes defined in this package itself
  (`MultiAssaySpatialExperiment`, `RasterLayerList`, `PointsLayerList`,
  `ShapesLayerList`, `SpatialLayerList`) were left as plain `\linkS4class{}`,
  since those never needed anchoring.
- Fixed the stale `http://` Bioconductor badge URL in `README.md` (`R CMD
  check` flagged it as a 301 redirect to `https://`).

## Bug fixes

- Added the missing `Matrix` `Suggests` dependency, used directly in
  `tests/` but not declared, which `R CMD check` flagged as an unstated
  test dependency.

# MultiAssaySpatialExperiment 0.99.9

## Documentation

- Fixed a real terminology collision the 0.99.7 copy edit missed: three
  vignettes built a shapes layer literally named `"regions"`, which collides
  with `spatialMap`'s own `region` column (the points/shapes layer name each
  row references). Renamed to `cells` in the introduction vignette and
  `anatomical_regions` in *use cases* (where the layer is tissue-scale, not
  per-cell); added a paragraph after `annotateWithRegions()` in the
  introduction distinguishing the pre-existing `region` column from the new
  annotation column the function adds. A `colData` column separately named
  `region` in *Working with MultiAssaySpatialExperiment*'s subsetting example
  collided for an unrelated reason and is now `zone`. The cheatsheet's
  `annotateWithRegions()` entry made the identical mistake in prose
  ("records the region on `spatialMap`") and now describes what the function
  actually does: add a new column named for the shapes layer.
- Introduced `ExperimentList`, `colData` and `sampleMap` before their first
  use in the introduction vignette instead of leaving them to appear
  unexplained in code output, added the "each row reads ..." sentence for
  `sampleMap` that *Design of MultiAssaySpatialExperiment* already used for
  `spatialMap`, explained `st_polygon()`/`st_sfc()`/`sfc_POLYGON` where they
  are first used, and noted that `read*()` functions exist for real platform
  data before walking through manual construction.
- Copy edit across all five vignettes for correctness, consistency, clarity
  and concision. Correctness: a cross-reference in *use cases* promised
  "coercion details and limitations" in a section that does not discuss
  them (removed); a heading and a summary bullet called a direct
  point-in-polygon count "the deconvolution result"/"cell type
  deconvolution" when the adjacent body text explicitly says it is not
  deconvolution (reworded to "composition result"/"cell type transfer");
  the *use cases* summary referenced a "deliberately plain R" sentence that
  had been removed from Workflow 2 earlier in this same round of edits
  (reworded to stand on its own); "FOV" was used before being spelled out;
  a comma splice in the `subsetByAssay()` note in *Working with
  MultiAssaySpatialExperiment*. Consistency: "Part 1: Data Import" and
  "Part 2: Real-World Workflows" used Title Case where every other heading
  in the same document uses sentence case after a colon; the Xenium reader
  section had an "Options" table that the Visium, CosMx and MERSCOPE
  sections did not (added, values taken from each function's own roxygen
  `@param` docs); two sentences in *use cases* used first-person-plural
  narrative voice ("We will create...") where the rest of the document
  describes what the code does. Concision: the entity-relationship diagram's
  intro sentence in *Design of MultiAssaySpatialExperiment* duplicated its
  own figure caption almost verbatim (trimmed); the "Scientific question"
  callout in *use cases* used more promotional language ("leverage...
  improve our understanding of") than the plain register used elsewhere in
  the document.

## Bug fixes

- `aggregateByRegion(x, by, FUN = "count")` named its result column `region`
  regardless of `by`; the `"sum"`/`"mean"` branches already name result
  columns after `by`, so `"count"` was the odd one out and, for any `by`
  other than a layer literally named `"region"`, produced a column name that
  did not match the argument that produced it. It now names the column after
  `by`, e.g. `by = "cells"` produces a `cells` column.

## Internal changes

- Replaced fourteen further `for` loops with `lapply()`/`vapply()`/
  `Reduce()`/vectorized code, in the same spirit as 0.99.8's refactor:
  `spatialMatch()`'s two loops (one needed no loop or `apply` at all, just
  vectorized arithmetic); `.subsetSpatialMapByColumn()` and
  `.subsetMASEBySpatialFilter()`; `.aggregateMatrixByRegion()` and
  `aggregateByRegion()`'s per-assay loop; `.mergeLayerList()` and
  `.reduceLayerList()` (the latter was a manual fold, now a `Reduce()`);
  four loops in `.check_visiumhd_files()`/`.build_mase_from_visiumhd()`, two
  of which re-tested a loop-invariant `is.null(boundaries)` condition on
  every iteration and are now hoisted out; the `ShapesLayerList()`
  constructor's coercion loop; the `annotGeometryNames` loop in the
  `SpatialFeatureExperiment` coercion; `.spatialImages_from_imgData()`'s
  rename-and-filter loop; and a test helper's loop around `file.remove()`,
  which already vectorizes. Verified against the full test suite plus
  targeted re-runs of every touched test file, and (for `spatialMatch()`'s
  vectorized branch) a direct comparison against the original loop's output.
- Consolidated the fifteen file-pattern-search loops kept as-is in 0.99.8
  (correctly: they `break` on the first match, so an apply variant can't
  short-circuit) into two shared helpers, `.firstExistingPath()` and
  `.firstMatchingFile()`, in a new "File discovery" section of
  `MultiAssaySpatialExperiment-read-helpers.R`. Those fifteen loops were the
  same six-line block copy-pasted across the CosMx, MERSCOPE, Xenium and
  Visium readers; all fifteen call sites now delegate to one of the two
  helpers instead, including the one Xenium case that also checks
  `dir.exists()` (passed in as a `check` override). Behavior is unchanged,
  confirmed by directly unit-testing both helpers against no-match, match,
  and the Xenium dual-check case, and by the full read-technology test file
  (151 assertions) passing unchanged. `for` loops in `R/` go from 43 to 15
  overall between this entry and 0.99.8; the remainder are the validity-check
  accumulators, side-effect-only loops, and upstream-fidelity ports (e.g.
  `.get_img_idx_mase()`) that 0.99.8 already explains.

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
