### =========================================================================
### MultiAssaySpatialExperiment - combining
### -------------------------------------------------------------------------
###
### c, cbind, complete.cases.
###
### -------------------------------------------------------------------------

#' @include MultiAssaySpatialExperiment-class.R
NULL

#' MultiAssaySpatialExperiment combining
#'
#' @description
#' Methods for combining \linkS4class{MultiAssaySpatialExperiment} objects.
#' \code{c} merges two MASE objects (adding assays); \code{cbind} combines along
#' columns (adding specimens); \code{complete.cases} identifies specimens with
#' data in all assays.
#'
#' @section Combining:
#' \describe{
#'   \item{\code{c(x, ..., sampleMap = NULL, mapFrom = NULL)}:}{
#'     When \code{...} is a single MultiAssaySpatialExperiment, merges the two
#'     objects (experiments, colData, sampleMap, spatial layers, imgData,
#'     spatialMap). Experiment and spatial layer names must be unique across
#'     objects. Otherwise delegates to \code{\link[MultiAssayExperiment]{c,MultiAssayExperiment-method}}.
#'   }
#'   \item{\code{cbind(..., deparse.level = 1)}:}{
#'     Combines MASE objects column-wise (same assay names, adds specimens).
#'     Spatial layers with the same name are row-bound; \code{imgData} and
#'     \code{spatialMap} are row-bound.
#'   }
#'   \item{\code{complete.cases(x)}:}{
#'     Returns a logical vector indicating which specimens (primaries) have
#'     data in all assays. Delegates to
#'     \code{\link[MultiAssayExperiment]{complete.cases,MultiAssayExperiment-method}}.
#'   }
#' }
#'
#' @author Patrick Aboyoun
#'
#' @aliases c,MultiAssaySpatialExperiment-method
#' @aliases cbind,MultiAssaySpatialExperiment-method
#'
#' @seealso
#' \itemize{
#'   \item \linkS4class{MultiAssaySpatialExperiment} for the main class
#'   \item \link{MultiAssaySpatialExperiment-subset} for subsetting methods
#' }
#'
#' @name MultiAssaySpatialExperiment-combine
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### c concatenation
###

#' @importFrom BiocGenerics cbind rbind
#' @importFrom MultiAssayExperiment experiments sampleMap ExperimentList
#' @importFrom S4Vectors DataFrame metadata merge
#' @importFrom SpatialExperiment imgData
#' @importFrom SummarizedExperiment colData assay
.mergeMASE <- function(x, y) {
    if (!is(y, "MultiAssaySpatialExperiment"))
        stop("'c' with two MultiAssaySpatialExperiment requires both to be ",
             "MultiAssaySpatialExperiment")
    if (any(names(x) %in% names(y)))
        stop("Provide unique experiment names")
    for (acc in list(spatialImages, spatialLabels, spatialPoints, spatialShapes)) {
        nms_x <- names(acc(x))
        nms_y <- names(acc(y))
        if (any(nms_x %in% nms_y))
            stop("Provide unique spatial layer names across objects")
    }
    expz <- c(experiments(x), experiments(y))
    sampz <- rbind(sampleMap(x), sampleMap(y))
    coldx <- colData(x)
    coldy <- colData(y)
    shared <- intersect(names(coldx), names(coldy))
    by <- if (length(shared)) c("row.names", shared) else "row.names"
    cdatz <- merge(x = coldx, y = coldy, by = by,
        all = TRUE, sort = FALSE, stringsAsFactors = FALSE)
    rownames(cdatz) <- cdatz[["Row.names"]]
    cdatz <- cdatz[, names(cdatz) != "Row.names", drop = FALSE]
    metaz <- c(metadata(x), metadata(y))
    imgx <- imgData(x)
    imgy <- imgData(y)
    new_img <- if (is.null(imgx) && is.null(imgy)) NULL
    else if (is.null(imgx)) imgy
    else if (is.null(imgy)) imgx
    else rbind(imgx, imgy)
    spx <- spatialMap(x)
    spy <- spatialMap(y)
    new_sp <- if (is.null(spx) && is.null(spy)) NULL
    else if (is.null(spx)) spy
    else if (is.null(spy)) spx
    else rbind(spx, spy)
    .mergeSpatialLayers <- function(acc) {
        lx <- acc(x)
        ly <- acc(y)
        if (length(lx) == 0L && length(ly) == 0L) return(acc(x)) # empty default
        nms <- c(names(lx), names(ly))
        lst <- c(as.list(lx), as.list(ly))
        cls <- class(lx)
        if (cls == "RasterLayerList")
            RasterLayerList(lst)
        else if (cls == "PointsLayerList")
            PointsLayerList(lst)
        else if (cls == "ShapesLayerList")
            ShapesLayerList(lst)
        else
            structure(lst, names = nms, class = class(lx))
    }
    imgs <- .mergeSpatialLayers(spatialImages)
    lbls <- .mergeSpatialLayers(spatialLabels)
    pts <- .mergeSpatialLayers(spatialPoints)
    shps <- .mergeSpatialLayers(spatialShapes)
    new("MultiAssaySpatialExperiment",
        ExperimentList = expz,
        colData = as(cdatz, "DataFrame"),
        sampleMap = sampz,
        metadata = metaz,
        images = imgs,
        labels = lbls,
        points = pts,
        shapes = shps,
        imgData = new_img,
        spatialMap = new_sp)
}

#' @exportMethod c
setMethod("c", "MultiAssaySpatialExperiment", function(x, ..., sampleMap = NULL, mapFrom = NULL) {
    args <- list(...)
    if (!length(args))
        stop("Provide experiments or a MultiAssaySpatialExperiment to concatenate")
    if (identical(length(args), 1L)) {
        input <- args[[1L]]
        if (is(input, "MultiAssaySpatialExperiment"))
            return(.mergeMASE(x, input))
    }
    callNextMethod()
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### cbind
###

.mergeLayerList <- function(lx, ly, constructor) {
    if (length(lx) == 0L && length(ly) == 0L)
        return(lx)
    if (length(lx) == 0L)
        return(ly)
    if (length(ly) == 0L)
        return(lx)
    all_nms <- union(names(lx), names(ly))
    lst <- setNames(vector("list", length(all_nms)), all_nms)
    for (nm in all_nms) {
        if (nm %in% names(lx) && nm %in% names(ly))
            lst[[nm]] <- rbind(lx[[nm]], ly[[nm]])
        else if (nm %in% names(lx))
            lst[[nm]] <- lx[[nm]]
        else
            lst[[nm]] <- ly[[nm]]
    }
    do.call(constructor, lst)
}

#' @exportMethod cbind
setMethod("cbind", "MultiAssaySpatialExperiment", function(..., deparse.level = 1) {
    args <- list(...)
    if (!length(args))
        stop("Provide at least one MultiAssaySpatialExperiment to cbind")
    if (!all(vapply(args, is, class2 = "MultiAssaySpatialExperiment", FUN.VALUE = TRUE)))
        stop("All arguments must be MultiAssaySpatialExperiment")
    if (length(args) == 1L)
        return(args[[1L]])
    nms <- unique(lapply(args, names))
    if (length(nms) > 1L)
        stop("All objects must have the same assay names")
    assay_nms <- nms[[1L]]
    if (length(assay_nms) == 0L)
        return(args[[1L]])
    new_exps <- ExperimentList(setNames(
        lapply(assay_nms, function(a) do.call(cbind, lapply(args, `[[`, a))),
        assay_nms))
    new_cd <- do.call(rbind, lapply(args, colData))
    new_sm <- do.call(rbind, lapply(args, sampleMap))
    new_meta <- do.call(c, unname(lapply(args, metadata)))
    img_list <- lapply(args, imgData)
    new_img <- if (all(vapply(img_list, is.null, logical(1L)))) NULL
    else do.call(rbind, Filter(Negate(is.null), img_list))
    sp_list <- lapply(args, spatialMap)
    new_sp <- if (all(vapply(sp_list, is.null, logical(1L)))) NULL
    else do.call(rbind, Filter(Negate(is.null), sp_list))
    .reduceLayerList <- function(acc, constructor) {
        layers <- lapply(args, acc)
        out <- layers[[1L]]
        for (i in seq(2L, length(layers)))
            out <- .mergeLayerList(out, layers[[i]], constructor)
        out
    }
    new_img_layers <- .reduceLayerList(spatialImages, RasterLayerList)
    new_lbl_layers <- .reduceLayerList(spatialLabels, RasterLayerList)
    new_pt_layers <- .reduceLayerList(spatialPoints, PointsLayerList)
    new_shp_layers <- .reduceLayerList(spatialShapes, ShapesLayerList)
    new("MultiAssaySpatialExperiment",
        ExperimentList = new_exps,
        colData = new_cd,
        sampleMap = new_sm,
        metadata = new_meta,
        images = new_img_layers,
        labels = new_lbl_layers,
        points = new_pt_layers,
        shapes = new_shp_layers,
        imgData = new_img,
        spatialMap = new_sp)
})

# complete.cases inherits from MultiAssayExperiment
