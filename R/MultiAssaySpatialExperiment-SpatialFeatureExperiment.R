### =========================================================================
### MultiAssaySpatialExperiment - SpatialFeatureExperiment integration
### -------------------------------------------------------------------------
###
### Coercion to SpatialFeatureExperiment; methods for dimGeometry, annotGeometry,
### spatialGraphs, bbox, unit, localResults, splitByCol, aggregate, findDebrisCells.
### Registered in .onLoad when SpatialFeatureExperiment is available.
###
### -------------------------------------------------------------------------

#' @include MultiAssaySpatialExperiment-class.R
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion (registered in .onLoad when SFE available)
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @importClassesFrom SpatialExperiment SpatialExperiment
#' @importFrom MultiAssayExperiment experiments
#' @importFrom SpatialExperiment imgData imgData<-
.register_SFE_coercion <- function() {
    setAs("MultiAssaySpatialExperiment", "SpatialFeatureExperiment",
        function(from) {
            exps <- experiments(from)
            if (length(exps) != 1L)
                stop("Coercion to SpatialFeatureExperiment requires exactly ",
                     "one assay; found ", length(exps))
            exp <- exps[[1L]]
            if (is(exp, "SpatialFeatureExperiment")) {
                sfe <- exp
            } else if (is(exp, "SpatialExperiment")) {
                sfe <- as(exp, "SpatialFeatureExperiment")
            } else {
                stop("The single assay must be a SpatialFeatureExperiment or ",
                     "SpatialExperiment; found ", class(exp)[1L])
            }
            mase_img <- imgData(from)
            if (!is.null(mase_img))
                imgData(sfe) <- mase_img
            sfe
        }
    )
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper: get SFE experiment for delegation
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @importClassesFrom SpatialExperiment SpatialExperiment
#' @importFrom MultiAssayExperiment experiments
.get_SFE_experiment <- function(x, assay = 1L) {
    exps <- experiments(x)
    if (length(exps) == 0L)
        return(NULL)
    exp <- exps[[assay]]
    if (is(exp, "SpatialFeatureExperiment"))
        exp
    else if (is(exp, "SpatialExperiment"))
        as(exp, "SpatialFeatureExperiment")
    else
        NULL
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helpers for splitByCol and aggregate
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @importFrom MultiAssayExperiment experiments
.sfe_to_mase <- function(sfe, template) {
    if (length(experiments(template)) != 1L)
        return(sfe)
    as(sfe, "MultiAssaySpatialExperiment")
}

#' @importFrom MultiAssayExperiment experiments
.sfe_list_to_mase_list <- function(sfe_list, template) {
    if (length(experiments(template)) != 1L)
        return(sfe_list)
    lapply(sfe_list, as, "MultiAssaySpatialExperiment")
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Methods (registered in .onLoad when SFE available)
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.register_SFE_methods <- function() {
    setMethod("dimGeometry", "MultiAssaySpatialExperiment",
        function(x, type = 1L, MARGIN, sample_id = 1L, withDimnames = TRUE) {
            exp <- .get_SFE_experiment(x)
            if (is.null(exp))
                stop("No SpatialFeatureExperiment assay in 'x'; ",
                     "geometry methods require an SFE experiment")
            SpatialFeatureExperiment::dimGeometry(exp,
                type = type, MARGIN = MARGIN, sample_id = sample_id,
                withDimnames = withDimnames)
        }
    )

    setMethod("dimGeometries", "MultiAssaySpatialExperiment",
        function(x, MARGIN = 2, withDimnames = TRUE) {
            exp <- .get_SFE_experiment(x)
            if (is.null(exp))
                stop("No SpatialFeatureExperiment assay in 'x'; ",
                     "geometry methods require an SFE experiment")
            SpatialFeatureExperiment::dimGeometries(exp,
                MARGIN = MARGIN, withDimnames = withDimnames)
        }
    )

    setMethod("dimGeometryNames", "MultiAssaySpatialExperiment",
        function(x, MARGIN) {
            exp <- .get_SFE_experiment(x)
            if (is.null(exp))
                stop("No SpatialFeatureExperiment assay in 'x'; ",
                     "geometry methods require an SFE experiment")
            SpatialFeatureExperiment::dimGeometryNames(exp, MARGIN)
        }
    )

    setMethod("annotGeometry", "MultiAssaySpatialExperiment",
        function(x, type = 1L, sample_id = 1L) {
            exp <- .get_SFE_experiment(x)
            if (is.null(exp))
                stop("No SpatialFeatureExperiment assay in 'x'; ",
                     "geometry methods require an SFE experiment")
            SpatialFeatureExperiment::annotGeometry(exp,
                type = type, sample_id = sample_id)
        }
    )

    setMethod("annotGeometries", "MultiAssaySpatialExperiment",
        function(x) {
            exp <- .get_SFE_experiment(x)
            if (is.null(exp))
                stop("No SpatialFeatureExperiment assay in 'x'; ",
                     "geometry methods require an SFE experiment")
            SpatialFeatureExperiment::annotGeometries(exp)
        }
    )

    setMethod("annotGeometryNames", "MultiAssaySpatialExperiment",
        function(x) {
            exp <- .get_SFE_experiment(x)
            if (is.null(exp))
                stop("No SpatialFeatureExperiment assay in 'x'; ",
                     "geometry methods require an SFE experiment")
            SpatialFeatureExperiment::annotGeometryNames(exp)
        }
    )

    setMethod("spatialGraphs", "MultiAssaySpatialExperiment",
        function(x, MARGIN = NULL, sample_id = "all", name = "all") {
            exp <- .get_SFE_experiment(x)
            if (is.null(exp))
                stop("No SpatialFeatureExperiment assay in 'x'; ",
                     "spatialGraphs requires an SFE experiment")
            SpatialFeatureExperiment::spatialGraphs(exp,
                MARGIN = MARGIN, sample_id = sample_id, name = name)
        }
    )

    setMethod("spatialGraph", "MultiAssaySpatialExperiment",
        function(x, type = 1L, MARGIN, sample_id = 1L) {
            exp <- .get_SFE_experiment(x)
            if (is.null(exp))
                stop("No SpatialFeatureExperiment assay in 'x'; ",
                     "spatialGraph requires an SFE experiment")
            SpatialFeatureExperiment::spatialGraph(exp,
                type = type, MARGIN = MARGIN, sample_id = sample_id)
        }
    )

    setMethod("spatialGraphNames", "MultiAssaySpatialExperiment",
        function(x, MARGIN, sample_id = 1L) {
            exp <- .get_SFE_experiment(x)
            if (is.null(exp))
                stop("No SpatialFeatureExperiment assay in 'x'; ",
                     "spatialGraphNames requires an SFE experiment")
            SpatialFeatureExperiment::spatialGraphNames(exp,
                MARGIN = MARGIN, sample_id = sample_id)
        }
    )

    setMethod("bbox", "MultiAssaySpatialExperiment",
        function(sfe, sample_id = "all", ...) {
            exp <- .get_SFE_experiment(sfe)
            if (is.null(exp))
                stop("No SpatialFeatureExperiment assay in 'x'; ",
                     "bbox requires an SFE experiment")
            SpatialFeatureExperiment::bbox(exp, sample_id = sample_id, ...)
        }
    )

    setMethod("unit", "MultiAssaySpatialExperiment",
        function(x) {
            exp <- .get_SFE_experiment(x)
            if (is.null(exp))
                return(NULL)
            SpatialFeatureExperiment::unit(exp)
        }
    )

    setMethod("localResults", "MultiAssaySpatialExperiment",
        function(x, sample_id = "all", name = "all", features = NULL,
                 colGeometryName = NULL, annotGeometryName = NULL,
                 withDimnames = TRUE, swap_rownames = NULL, ...) {
            exp <- .get_SFE_experiment(x)
            if (is.null(exp))
                stop("No SpatialFeatureExperiment assay in 'x'; ",
                     "localResults requires an SFE experiment")
            SpatialFeatureExperiment::localResults(exp,
                sample_id = sample_id, name = name, features = features,
                colGeometryName = colGeometryName,
                annotGeometryName = annotGeometryName,
                withDimnames = withDimnames,
                swap_rownames = swap_rownames, ...)
        }
    )

    setMethod("localResultNames", "MultiAssaySpatialExperiment",
        function(x) {
            exp <- .get_SFE_experiment(x)
            if (is.null(exp))
                return(character(0L))
            SpatialFeatureExperiment::localResultNames(exp)
        }
    )

    setMethod("splitByCol", c("MultiAssaySpatialExperiment", "sf"),
        function(x, f, sample_id = "all", colGeometryName = 1L, cover = FALSE) {
            exp <- .get_SFE_experiment(x)
            if (is.null(exp))
                stop("No SpatialFeatureExperiment assay in 'x'; ",
                     "splitByCol requires an SFE experiment")
            spl <- SpatialFeatureExperiment::splitByCol(exp,
                f = f, sample_id = sample_id,
                colGeometryName = colGeometryName, cover = cover)
            .sfe_list_to_mase_list(spl, x)
        }
    )

    setMethod("splitByCol", c("MultiAssaySpatialExperiment", "sfc"),
        function(x, f, sample_id = 1L, colGeometryName = 1L, cover = FALSE) {
            exp <- .get_SFE_experiment(x)
            if (is.null(exp))
                stop("No SpatialFeatureExperiment assay in 'x'; ",
                     "splitByCol requires an SFE experiment")
            spl <- SpatialFeatureExperiment::splitByCol(exp,
                f = f, sample_id = sample_id,
                colGeometryName = colGeometryName, cover = cover)
            .sfe_list_to_mase_list(spl, x)
        }
    )

    setMethod("splitByCol", c("MultiAssaySpatialExperiment", "list"),
        function(x, f, sample_id = "all", colGeometryName = 1L, cover = FALSE) {
            exp <- .get_SFE_experiment(x)
            if (is.null(exp))
                stop("No SpatialFeatureExperiment assay in 'x'; ",
                     "splitByCol requires an SFE experiment")
            spl <- SpatialFeatureExperiment::splitByCol(exp,
                f = f, sample_id = sample_id,
                colGeometryName = colGeometryName, cover = cover)
            .sfe_list_to_mase_list(spl, x)
        }
    )

    setMethod("aggregate", "MultiAssaySpatialExperiment",
        function(x, by = NULL, sample_id = "all", ...) {
            exp <- .get_SFE_experiment(x)
            if (is.null(exp))
                stop("No SpatialFeatureExperiment assay in 'x'; ",
                     "aggregate requires an SFE experiment")
            agg <- SpatialFeatureExperiment::aggregate(exp,
                by = by, sample_id = sample_id, ...)
            .sfe_to_mase(agg, x)
        }
    )

    setMethod("findDebrisCells", "MultiAssaySpatialExperiment",
        function(x, max_cells = 5, distance_cutoff = 50,
                 BNPARAM = NULL, BPPARAM = BiocParallel::SerialParam()) {
            exp <- .get_SFE_experiment(x)
            if (is.null(exp))
                stop("No SpatialFeatureExperiment assay in 'x'; ",
                     "findDebrisCells requires an SFE experiment")
            SpatialFeatureExperiment::findDebrisCells(exp,
                max_cells = max_cells, distance_cutoff = distance_cutoff,
                BNPARAM = BNPARAM, BPPARAM = BPPARAM)
        }
    )
}
