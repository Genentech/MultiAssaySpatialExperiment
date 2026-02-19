### =========================================================================
### aggregateByRegion — multi-assay spatial aggregation
### -------------------------------------------------------------------------
###
### Aggregate assay values by shape region. Requires annotateWithRegions first.
###
### -------------------------------------------------------------------------

#' @include MultiAssaySpatialExperiment-class.R
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### aggregateByRegion
###

#' Aggregate assay values by shape region
#'
#' @description
#' Aggregate assay values (e.g., counts, expression) by shape region. Each assay
#' column is linked to a point via \code{spatialMap}; the annotation column
#' (from \code{\link{annotateWithRegions}}) links points to shapes. Values are
#' aggregated within each shape.
#'
#' @param x A \linkS4class{MultiAssaySpatialExperiment}.
#' @param by Character. Name of the region column in \code{spatialMap} that holds
#'     the shape instance IDs (the \code{regionCol} from
#'     \code{\link{annotateWithRegions}}, default is the shapes layer name).
#' @param assays Character. Assay names to aggregate (default: all assays in
#'     \code{spatialMap}).
#' @param FUN Character. Aggregation function: \code{"sum"}, \code{"mean"},
#'     or \code{"count"}. For \code{"count"}, returns the number of points per
#'     shape; for \code{"sum"} and \code{"mean"}, aggregates assay values.
#'
#' @details
#' Run \code{\link{annotateWithRegions}} first to add the \code{by} column to
#' \code{spatialMap}. Rows with \code{NA} in \code{by} are dropped.
#'
#' For \code{"sum"} and \code{"mean"}, the assay matrix (features x columns) is
#' grouped by shape: columns that map to the same shape are aggregated. The
#' result is a list of matrices (one per assay), each with features as rows and
#' shape instance IDs as columns.
#'
#' @return
#' A list of matrices (or for \code{FUN = "count"} a single \code{DataFrame}),
#' one element per assay. Each matrix has features (rows) and shape instance IDs
#' (columns).
#'
#' @section Polymorphism:
#' Uses S3/S4 generics only. Table-like elements (\code{spatialMap}, experiments)
#' use \code{nrow}, \code{colnames}, \code{rownames}, \code{[[}, \code{[},
#' \code{unique}, and \code{split}. Assay data use \code{\link[SummarizedExperiment]{assay}}
#' for SummarizedExperiment-like objects; plain matrices are passed through.
#' Layer lists use \code{names}. Consistent with
#' \code{\link{annotateWithRegions}} and \code{\link{subsetByPolygon}}; see the
#' Subset vignette.
#'
#' @author Patrick Aboyoun
#'
#' @examples
#' pts <- S4Vectors::DataFrame(
#'     x = c(1.5, 2.5, 3.5),
#'     y = c(1.5, 2.5, 3.5),
#'     instance_id = c("A", "B", "C"))
#' shp_df <- S4Vectors::DataFrame(
#'     instance_id = c("cell1", "cell2", "cell3"),
#'     geometry = sf::st_sfc(
#'         sf::st_polygon(list(matrix(c(1,1,2,1,2,2,1,2,1,1), ncol=2, byrow=TRUE))),
#'         sf::st_polygon(list(matrix(c(2,1,3,1,3,2,2,2,2,1), ncol=2, byrow=TRUE))),
#'         sf::st_polygon(list(matrix(c(2,2,3,2,3,3,2,3,2,2), ncol=2, byrow=TRUE)))))
#' mase <- MultiAssaySpatialExperiment(
#'     experiments = ExperimentList(assay1 = matrix(1:9, 3, 3,
#'         dimnames = list(paste0("G", 1:3), c("A", "B", "C")))),
#'     colData = S4Vectors::DataFrame(row.names = "s1"),
#'     sampleMap = S4Vectors::DataFrame(
#'         assay = "assay1", primary = "s1", colname = c("A", "B", "C")),
#'     points = PointsLayerList(centroids = pts),
#'     shapes = ShapesLayerList(cells = shp_df),
#'     spatialMap = S4Vectors::DataFrame(
#'         assay = "assay1", colname = c("A", "B", "C"),
#'         region = "centroids", instance_id = c("A", "B", "C")))
#' mase <- annotateWithRegions(mase, points = "centroids", shapes = "cells")
#' agg <- aggregateByRegion(mase, by = "cells", FUN = "sum")
#' agg[["assay1"]]
#'
#' @aliases
#' aggregateByRegion
#' aggregateByRegion,MultiAssaySpatialExperiment-method
#'
#' @export
setGeneric("aggregateByRegion",
    function(x, by, assays = NULL, FUN = c("sum", "mean", "count"))
        standardGeneric("aggregateByRegion"))

#' @importFrom SummarizedExperiment assay
.getAssayMatrix <- function(ex) {
    if (is.matrix(ex)) ex else assay(ex, 1L)
}

.aggregateMatrixByRegion <- function(mat, sm, by, FUN) {
    if (length(dim(mat)) != 2L)
        return(NULL)
    cn <- colnames(mat)
    if (is.null(cn))
        return(NULL)
    rn <- rownames(mat)
    if (is.null(rn))
        rn <- seq_len(nrow(mat))
    if (!"colname" %in% colnames(sm) || !by %in% colnames(sm))
        return(NULL)
    col_by_region <- split(sm[["colname"]], sm[[by]])
    regions <- names(col_by_region)
    out <- matrix(NA_real_, nrow = nrow(mat), ncol = length(regions),
        dimnames = list(rn, regions))
    for (reg in regions) {
        cols <- intersect(col_by_region[[reg]], cn)
        if (length(cols) == 0L)
            next
        sub <- mat[, cols, drop = FALSE]
        out[, reg] <- if (FUN == "sum") rowSums(sub) else rowMeans(sub)
    }
    out
}

#' @importFrom MultiAssayExperiment experiments
#' @importFrom S4Vectors DataFrame
#' @exportMethod aggregateByRegion
setMethod("aggregateByRegion", "MultiAssaySpatialExperiment",
function(x, by, assays = NULL, FUN = c("sum", "mean", "count")) {
    FUN <- match.arg(FUN)
    spmap <- spatialMap(x)
    if (is.null(spmap) || nrow(spmap) == 0L)
        return(list())
    if (!by %in% colnames(spmap))
        stop("'by' column '", by, "' not found in spatialMap; run annotateWithRegions first")
    valid <- !is.na(spmap[[by]])
    if (!any(valid))
        return(list())
    spmap <- spmap[valid, , drop = FALSE]
    assay_nms <- if (is.null(assays)) unique(spmap[["assay"]]) else assays
    assay_nms <- intersect(assay_nms, names(experiments(x)))
    if (length(assay_nms) == 0L)
        return(list())

    if (FUN == "count") {
        cnt <- table(spmap[[by]])
        return(DataFrame(region = names(cnt), count = as.integer(cnt)))
    }

    result <- list()
    exps <- experiments(x)
    for (nm in assay_nms) {
        sm <- spmap[spmap[["assay"]] == nm, , drop = FALSE]
        if (nrow(sm) == 0L)
            next
        ex <- exps[[nm]]
        mat <- .getAssayMatrix(ex)
        out <- .aggregateMatrixByRegion(mat, sm, by, FUN)
        if (!is.null(out))
            result[[nm]] <- out
    }
    result
})
