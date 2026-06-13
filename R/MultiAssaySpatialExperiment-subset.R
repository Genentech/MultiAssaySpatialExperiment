### =========================================================================
### MultiAssaySpatialExperiment - subsetting
### -------------------------------------------------------------------------
###
### subsetByColData, subsetByRow, subsetByColumn, subsetByAssay, subsetByBoundingBox,
### subsetByPolygon, and [ with propagation to spatial slots (spatialMap,
### imgData, images, labels, points, shapes).
###
### -------------------------------------------------------------------------

#' @include MultiAssaySpatialExperiment-class.R
NULL

#' MultiAssaySpatialExperiment subsetting
#'
#' @description
#' Methods for subsetting a \linkS4class{MultiAssaySpatialExperiment}. They
#' extend \code{\link[MultiAssayExperiment]{subsetBy}} with propagation to
#' spatial slots.
#'
#' @section Subsetting:
#' In the snippets below, \code{x} is a
#' \linkS4class{MultiAssaySpatialExperiment}:
#' \describe{
#'   \item{\code{subsetByColData(x, y)}:}{
#'     Subset by primary identifiers (specimens). \code{spatialMap} and
#'     \code{imgData} are filtered to retained specimens; point and shape rows
#'     are not trimmed unless removed via the map.
#'   }
#'   \item{\code{subsetByRow(x, y, i = TRUE, ...)}:}{
#'     Subset rows of experiments. Spatial layers are assay-level and unchanged.
#'   }
#'   \item{\code{subsetByColumn(x, y)}:}{
#'     Subset assay columns (observations). When \code{y} is a \code{list}
#'     (one element per assay), each element may be column names, indices, or
#'     a logical vector — the usual \code{\link[MultiAssayExperiment]{subsetByColumn}}
#'     interface. \code{spatialMap}, \code{imgData}, and linked points, shapes,
#'     images, and labels are filtered to retained columns. To subset by assay
#'     \code{colData}, build a logical vector or column names and pass a list,
#'     e.g. \code{subsetByColumn(mase, list(rna = colData(experiments(mase)[["rna"]])$tissue == "tumor"))}.
#'   }
#'   \item{\code{subsetByAssay(x, y)}:}{
#'     Subset assays. Points, shapes, images, labels, \code{spatialMap}, and
#'     \code{imgData} are filtered to layers and rows linked to retained assays
#'     via \code{spatialMap}. Auxiliary shape layers not referenced in
#'     \code{spatialMap} are dropped.
#'   }
#'   \item{\code{x[i, j, k, ..., drop = FALSE]}:}{
#'     Subset by row (\code{i}), column (\code{j}), and assay (\code{k}) indices.
#'     Equivalent to composing \code{subsetByColData}, \code{subsetByColumn},
#'     \code{subsetByRow}, and \code{subsetByAssay}. If \code{drop = TRUE},
#'     empty assays are removed.
#'   }
#'   \item{\code{subsetByBoundingBox(x, xmin, xmax, ymin, ymax, ...)}:}{
#'     Subset by bounding box. Filters points and shapes within the rectangle;
#'     propagates to assays via \code{spatialMap} when present.
#'   }
#'   \item{\code{subsetByPolygon(x, polygon, ...)}:}{
#'     Subset by polygon. Filters points and shapes within or intersecting the
#'     \code{sf} geometry; propagates to assays via \code{spatialMap} when present.
#'   }
#' }
#'
#' @section Polymorphism:
#' All subsetting uses S3/S4 generics only. Table-like elements
#' (\code{spatialMap}, \code{sampleMap}, \code{imgData}, points, shapes) use
#' \code{nrow}, \code{colnames}, \code{[[}, \code{[}, \code{unique}, and
#' \code{split} (from \pkg{BiocGenerics}). Layer lists use
#' \code{names} and \code{length}. Spatial filtering additionally requires
#' \code{sf::st_intersects} (and \code{st_as_sfc}, \code{st_bbox}
#' for shapes / polygon). For spatial filtering, \code{instance_id} in points,
#' shapes, and \code{spatialMap} must have matching types (no coercion). See
#' vignette \emph{Working with MultiAssaySpatialExperiment} for a per-operation
#' propagation table.
#'
#' @author Patrick Aboyoun
#'
#' @aliases subsetByColData,MultiAssaySpatialExperiment,ANY-method
#' @aliases subsetByColData,MultiAssaySpatialExperiment,character-method
#' @aliases subsetByRow,MultiAssaySpatialExperiment,ANY-method
#' @aliases subsetByRow,MultiAssaySpatialExperiment,list-method
#' @aliases subsetByRow,MultiAssaySpatialExperiment,List-method
#' @aliases subsetByColumn,MultiAssaySpatialExperiment,ANY-method
#' @aliases subsetByAssay,MultiAssaySpatialExperiment-method
#' @aliases subsetByBoundingBox
#' @aliases subsetByBoundingBox,MultiAssaySpatialExperiment-method
#' @aliases subsetByPolygon
#' @aliases subsetByPolygon,MultiAssaySpatialExperiment-method
#' @aliases [,MultiAssaySpatialExperiment,ANY,ANY,ANY-method
#'
#' @seealso
#' \itemize{
#'   \item \linkS4class{MultiAssaySpatialExperiment} for the main class
#'   \item \link{MultiAssaySpatialExperiment-combine} for combining methods
#' }
#'
#' @keywords methods
#'
#' @name MultiAssaySpatialExperiment-subset
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### subsetBy methods
###

.subsetSpatialMapByColData <- function(spatialMap, sampleMap) {
    if (is.null(spatialMap) || nrow(spatialMap) == 0L)
        return(spatialMap)
    cols <- c("assay", "colname")
    if (!all(cols %in% colnames(spatialMap)))
        return(spatialMap)
    sm_cols <- sampleMap[, c("assay", "colname"), drop = FALSE]
    sm_keys <- paste(sm_cols[["assay"]], sm_cols[["colname"]], sep = "\r")
    sp_keys <- paste(spatialMap[["assay"]], spatialMap[["colname"]], sep = "\r")
    spatialMap[sp_keys %in% sm_keys, , drop = FALSE]
}

.subsetImgDataByColData <- function(imgData, primary_ids) {
    if (is.null(imgData) || nrow(imgData) == 0L)
        return(imgData)
    if (!"sample_id" %in% colnames(imgData))
        return(imgData)
    imgData[imgData[["sample_id"]] %in% primary_ids, , drop = FALSE]
}

.subsetSpatialLayersByRegion <- function(elements, spatialMap, element_type) {
    nms <- names(elements)
    if (length(nms) == 0L)
        return(elements)
    if (is.null(spatialMap) || nrow(spatialMap) == 0L)
        return(elements)
    if (!"region" %in% colnames(spatialMap) || !"element_type" %in% colnames(spatialMap))
        return(elements)

    ## Filter spatialMap to rows for this element_type only
    et_rows <- spatialMap[["element_type"]] == element_type
    if (!any(et_rows))
        return(elements[integer(0L)])

    regions <- unique(spatialMap[["region"]][et_rows])
    keep <- nms %in% regions
    if (!any(keep))
        return(elements[integer(0L)])

    elements[keep]
}

.subsetSpatialLayersByAssay <- function(elements, assay_names) {
    nms <- names(elements)
    if (length(nms) == 0L)
        return(elements)
    keep <- nms %in% assay_names
    if (!any(keep))
        return(elements)
    elements[keep]
}

.subsetSpatialMapByAssay <- function(spatialMap, assay_names) {
    if (is.null(spatialMap) || nrow(spatialMap) == 0L)
        return(spatialMap)
    if (!"assay" %in% colnames(spatialMap))
        return(spatialMap)
    spatialMap[spatialMap[["assay"]] %in% assay_names, , drop = FALSE]
}

.subsetImgDataByAssay <- function(imgData, sampleMap, assay_names) {
    if (is.null(imgData) || nrow(imgData) == 0L)
        return(imgData)
    if (is.null(sampleMap) || nrow(sampleMap) == 0L)
        return(imgData)
    sm_sub <- sampleMap[sampleMap[["assay"]] %in% assay_names, , drop = FALSE]
    primary_ids <- unique(sm_sub[["primary"]])
    .subsetImgDataByColData(imgData, primary_ids)
}

.subsetSpatialMapByColumn <- function(spatialMap, experiments) {
    if (is.null(spatialMap) || nrow(spatialMap) == 0L)
        return(spatialMap)
    cols <- c("assay", "colname")
    if (!all(cols %in% colnames(spatialMap)))
        return(spatialMap)
    exp_names <- names(experiments)
    kept <- character()
    for (nm in exp_names) {
        cn <- colnames(experiments[[nm]])
        kept <- c(kept, paste(nm, cn, sep = "\r"))
    }
    sp_keys <- paste(spatialMap[["assay"]], spatialMap[["colname"]], sep = "\r")
    spatialMap[sp_keys %in% kept, , drop = FALSE]
}

.subsetSpatialElementsBySpatialMap <- function(x, spmap) {
    regions <- character()
    if (!is.null(spmap) && nrow(spmap) > 0L && "region" %in% colnames(spmap))
        regions <- unique(as.character(spmap[["region"]]))
    pts <- spatialPoints(x)
    shps <- spatialShapes(x)
    imgs <- spatialImages(x)
    labs <- spatialLabels(x)
    if (length(regions)) {
        drop_pts <- setdiff(names(pts), regions)
        if (length(drop_pts))
            pts <- pts[setdiff(seq_along(pts), match(drop_pts, names(pts)))]
        drop_shps <- setdiff(names(shps), regions)
        if (length(drop_shps))
            shps <- shps[setdiff(seq_along(shps), match(drop_shps, names(shps)))]
        drop_imgs <- setdiff(names(imgs), regions)
        if (length(drop_imgs))
            imgs <- imgs[setdiff(seq_along(imgs), match(drop_imgs, names(imgs)))]
        drop_labs <- setdiff(names(labs), regions)
        if (length(drop_labs))
            labs <- labs[setdiff(seq_along(labs), match(drop_labs, names(labs)))]
    }
    inst_ids <- if (!is.null(spmap) && nrow(spmap) > 0L &&
            "instance_id" %in% colnames(spmap))
        unique(spmap[["instance_id"]]) else character()
    for (reg in regions) {
        if (reg %in% names(pts)) {
            el <- pts[[reg]]
            if ("instance_id" %in% colnames(el) && length(inst_ids))
                pts[[reg]] <- el[as.character(el[["instance_id"]]) %in%
                    as.character(inst_ids), , drop = FALSE]
        }
        if (reg %in% names(shps)) {
            el <- shps[[reg]]
            if ("instance_id" %in% colnames(el) && length(inst_ids))
                shps[[reg]] <- el[as.character(el[["instance_id"]]) %in%
                    as.character(inst_ids), , drop = FALSE]
        }
    }
    replaceSlots(x,
        points = PointsLayerList(as.list(pts)),
        shapes = ShapesLayerList(as.list(shps)),
        images = RasterLayerList(as.list(imgs)),
        labels = RasterLayerList(as.list(labs)),
        check = FALSE)
}

.subsetPointsByBbox <- function(pt, xmin, xmax, ymin, ymax, x_col, y_col) {
    if (is.null(pt) || nrow(pt) == 0L)
        return(integer(0L))
    if (!x_col %in% colnames(pt) || !y_col %in% colnames(pt))
        return(integer(0L))
    xv <- pt[[x_col]]
    yv <- pt[[y_col]]
    keep <- (xv >= xmin & xv <= xmax) & (yv >= ymin & yv <= ymax)
    which(keep)
}

.subsetPointsByPolygon <- function(pt, polygon, x_col, y_col) {
    if (is.null(pt) || nrow(pt) == 0L)
        return(integer(0L))
    if (!x_col %in% colnames(pt) || !y_col %in% colnames(pt))
        return(integer(0L))
    which(spatialOverlaps(pt, polygon, coords = c(x_col, y_col)))
}

#' @importFrom sf st_as_sfc st_bbox
.subsetShapesByBbox <- function(shp, xmin, xmax, ymin, ymax) {
    if (is.null(shp) || nrow(shp) == 0L)
        return(integer(0L))
    if (!"geometry" %in% colnames(shp))
        return(integer(0L))
    bbox_sfc <- st_as_sfc(st_bbox(c(xmin = xmin, ymin = ymin,
                                     xmax = xmax, ymax = ymax)))
    which(spatialOverlaps(shp, bbox_sfc, geom = "geometry"))
}

.subsetShapesByPolygon <- function(shp, polygon) {
    if (is.null(shp) || nrow(shp) == 0L)
        return(integer(0L))
    if (!"geometry" %in% colnames(shp))
        return(integer(0L))
    which(spatialOverlaps(shp, polygon, geom = "geometry"))
}

.getInstanceIds <- function(el, idx) {
    if (length(idx) == 0L)
        return(character(0L))
    if ("instance_id" %in% colnames(el))
        el[["instance_id"]][idx]
    else
        idx
}

#' @importFrom MultiAssayExperiment experiments sampleMap
.subsetMASEBySpatialFilter <- function(x, retained_by_region) {
    if (length(retained_by_region) == 0L)
        return(x)
    spmap <- spatialMap(x)
    if (is.null(spmap) || nrow(spmap) == 0L) {
        warning("No spatialMap; assays unchanged; only points and shapes filtered",
                call. = FALSE)
        return(x)
    }
    if (!all(c("assay", "colname", "element_type", "region", "instance_id") %in% colnames(spmap)))
        return(x)
    keep_sp <- logical(nrow(spmap))
    for (i in seq_len(nrow(spmap))) {
        reg <- spmap[i, "region"]
        inst <- spmap[i, "instance_id"]
        kept_inst <- retained_by_region[[reg]]
        if (!is.null(kept_inst) && inst %in% kept_inst)
            keep_sp[[i]] <- TRUE
    }
    spmap <- spmap[keep_sp, , drop = FALSE]
    if (nrow(spmap) == 0L)
        return(x)
    j_list <- split(spmap[["colname"]], spmap[["assay"]])
    exps <- experiments(x)
    assay_nms <- names(exps)
    j_list_full <- lapply(assay_nms, function(a) {
        if (a %in% names(j_list))
            j_list[[a]]
        else
            colnames(exps[[a]])
    })
    names(j_list_full) <- assay_nms
    subsetByColumn(x, j_list_full)
}

#' @importFrom MultiAssayExperiment sampleMap subsetByColData
#' @importFrom SpatialExperiment imgData
#' @importFrom SummarizedExperiment colData
#' @exportMethod subsetByColData
setMethod("subsetByColData", c("MultiAssaySpatialExperiment", "ANY"),
function(x, y) {
    x <- callNextMethod()
    sm <- sampleMap(x)
    cd <- colData(x)
    replaceSlots(x,
                 spatialMap = .subsetSpatialMapByColData(spatialMap(x), sm),
                 imgData = .subsetImgDataByColData(imgData(x), rownames(cd)),
                 check = FALSE)
})

#' @importFrom MultiAssayExperiment subsetByColData sampleMap
#' @importFrom SummarizedExperiment colData
#' @importFrom SpatialExperiment imgData
#' @exportMethod subsetByColData
setMethod("subsetByColData", c("MultiAssaySpatialExperiment", "character"),
function(x, y) {
    x <- callNextMethod()
    sm <- sampleMap(x)
    cd <- colData(x)
    replaceSlots(x,
                 spatialMap = .subsetSpatialMapByColData(spatialMap(x), sm),
                 imgData = .subsetImgDataByColData(imgData(x), rownames(cd)),
                 check = FALSE)
})

#' @importFrom MultiAssayExperiment subsetByRow
#' @exportMethod subsetByRow
setMethod("subsetByRow", c("MultiAssaySpatialExperiment", "ANY"),
    function(x, y, i = TRUE, ...) callNextMethod())

#' @importFrom MultiAssayExperiment subsetByRow
#' @exportMethod subsetByRow
setMethod("subsetByRow", c("MultiAssaySpatialExperiment", "list"),
    function(x, y, ...) callNextMethod())

#' @importFrom MultiAssayExperiment subsetByRow
#' @exportMethod subsetByRow
setMethod("subsetByRow", c("MultiAssaySpatialExperiment", "List"),
    function(x, y, ...) callNextMethod())

#' @importFrom MultiAssayExperiment experiments sampleMap subsetByColumn
#' @importFrom SpatialExperiment imgData
#' @importFrom SummarizedExperiment colData
#' @exportMethod subsetByColumn
setMethod("subsetByColumn", c("MultiAssaySpatialExperiment", "ANY"),
function(x, y) {
    if (is.character(y) || is.logical(y) || is.numeric(y)) {
        callNextMethod()
    } else {
        x <- callNextMethod()
        cd <- colData(x)
        spmap <- .subsetSpatialMapByColumn(spatialMap(x), experiments(x))
        x <- replaceSlots(x,
            spatialMap = spmap,
            imgData = .subsetImgDataByColData(imgData(x), rownames(cd)),
            check = FALSE)
        .subsetSpatialElementsBySpatialMap(x, spmap)
    }
})

#' @importFrom MultiAssayExperiment experiments sampleMap subsetByAssay
#' @importFrom SpatialExperiment imgData
#' @exportMethod subsetByAssay
setMethod("subsetByAssay", c("MultiAssaySpatialExperiment", "ANY"),
function(x, y) {
    x <- callNextMethod()
    assay_names <- names(experiments(x))
    spmap <- .subsetSpatialMapByAssay(spatialMap(x), assay_names)
    replaceSlots(x,
                 images = .subsetSpatialLayersByAssay(spatialImages(x),
                                                     assay_names),
                 labels = .subsetSpatialLayersByAssay(spatialLabels(x),
                                                      assay_names),
                 points = .subsetSpatialLayersByRegion(spatialPoints(x), spmap, "points"),
                 shapes = .subsetSpatialLayersByRegion(spatialShapes(x), spmap, "shapes"),
                 spatialMap = spmap,
                 imgData = .subsetImgDataByAssay(imgData(x), sampleMap(x),
                                                 assay_names),
                 check = FALSE)
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### subsetByBoundingBox / subsetByPolygon
###

#' @export
setGeneric("subsetByBoundingBox", function(x, xmin, xmax, ymin, ymax, ...)
    standardGeneric("subsetByBoundingBox"))

#' @export
setGeneric("subsetByPolygon", function(x, polygon, ...)
    standardGeneric("subsetByPolygon"))

#' @importFrom MultiAssayExperiment experiments
#' @importFrom S4Vectors DataFrame
#' @exportMethod subsetByBoundingBox
setMethod("subsetByBoundingBox", "MultiAssaySpatialExperiment",
function(x, xmin, xmax, ymin, ymax, x_col = "x", y_col = "y", ...) {
    pts <- spatialPoints(x)
    shps <- spatialShapes(x)
    retained_by_region <- list()
    pt_nms <- names(pts)
    for (nm in pt_nms) {
        el <- pts[[nm]]
        idx <- .subsetPointsByBbox(el, xmin, xmax, ymin, ymax, x_col, y_col)
        retained_by_region[[nm]] <- .getInstanceIds(el, idx)
    }
    shp_nms <- names(shps)
    for (nm in shp_nms) {
        el <- shps[[nm]]
        idx <- .subsetShapesByBbox(el, xmin, xmax, ymin, ymax)
        if (length(idx) > 0L)
            retained_by_region[[nm]] <- .getInstanceIds(el, idx)
    }
    x <- .subsetMASEBySpatialFilter(x, retained_by_region)
    filt_pts <- list()
    for (nm in pt_nms) {
        el <- pts[[nm]]
        idx <- .subsetPointsByBbox(el, xmin, xmax, ymin, ymax, x_col, y_col)
        if (length(idx) > 0L)
            filt_pts[[nm]] <- el[idx, , drop = FALSE]
    }
    filt_shps <- list()
    for (nm in shp_nms) {
        el <- shps[[nm]]
        idx <- .subsetShapesByBbox(el, xmin, xmax, ymin, ymax)
        if (length(idx) > 0L)
            filt_shps[[nm]] <- el[idx, , drop = FALSE]
    }
    replaceSlots(x,
                points = PointsLayerList(filt_pts),
                shapes = ShapesLayerList(filt_shps),
                check = FALSE)
})

#' @importFrom MultiAssayExperiment experiments
#' @importFrom S4Vectors DataFrame
#' @importFrom sf st_sfc st_geometry
#' @exportMethod subsetByPolygon
setMethod("subsetByPolygon", "MultiAssaySpatialExperiment",
function(x, polygon, x_col = "x", y_col = "y", ...) {
    if (inherits(polygon, "sfg"))
        polygon <- st_sfc(polygon)
    else if (inherits(polygon, "sf"))
        polygon <- st_geometry(polygon)
    pts <- spatialPoints(x)
    shps <- spatialShapes(x)
    retained_by_region <- list()
    pt_nms <- names(pts)
    for (nm in pt_nms) {
        el <- pts[[nm]]
        idx <- .subsetPointsByPolygon(el, polygon, x_col, y_col)
        retained_by_region[[nm]] <- .getInstanceIds(el, idx)
    }
    shp_nms <- names(shps)
    for (nm in shp_nms) {
        el <- shps[[nm]]
        idx <- .subsetShapesByPolygon(el, polygon)
        if (length(idx) > 0L)
            retained_by_region[[nm]] <- .getInstanceIds(el, idx)
    }
    x <- .subsetMASEBySpatialFilter(x, retained_by_region)
    filt_pts <- list()
    for (nm in pt_nms) {
        el <- pts[[nm]]
        idx <- .subsetPointsByPolygon(el, polygon, x_col, y_col)
        if (length(idx) > 0L)
            filt_pts[[nm]] <- el[idx, , drop = FALSE]
    }
    filt_shps <- list()
    for (nm in shp_nms) {
        el <- shps[[nm]]
        idx <- .subsetShapesByPolygon(el, polygon)
        if (length(idx) > 0L)
            filt_shps[[nm]] <- el[idx, , drop = FALSE]
    }
    replaceSlots(x,
                points = PointsLayerList(filt_pts),
                shapes = ShapesLayerList(filt_shps),
                check = FALSE)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###

#' @importFrom MultiAssayExperiment ExperimentList experiments experiments<-
#' @importFrom MultiAssayExperiment drops<- subsetByAssay
.dropEmptyMASE <- function(object, warn = TRUE) {
    exps <- experiments(object)
    isEmptyAssay <- vapply(exps, function(e) {
        d <- dim(e)
        isTRUE((d[1L] == 0L) || (d[2L] == 0L))
    }, logical(1L))
    if (all(isEmptyAssay)) {
        drops(object) <- list(experiments = names(object))
        if (warn)
            warning("'experiments' dropped; see 'drops()'", call. = FALSE)
        experiments(object) <- ExperimentList()
    } else if (any(isEmptyAssay)) {
        keeps <- names(isEmptyAssay)[!isEmptyAssay]
        drops(object) <- list(experiments = names(isEmptyAssay)[isEmptyAssay])
        if (warn)
            warning("'experiments' dropped; see 'drops()'", call. = FALSE)
        object <- subsetByAssay(object, keeps)
    }
    object
}

#' @importFrom MultiAssayExperiment subsetByAssay subsetByColData subsetByColumn
#' @importFrom MultiAssayExperiment subsetByRow
.subsetMultiAssaySpatialExperiment <- function(x, i, j, k, ..., drop = FALSE) {
    if (missing(i) && missing(j) && missing(k)) {
        return(x)
    }
    if (!missing(j)) {
        if (is(j, "list") || is(j, "List"))
            x <- subsetByColumn(x, j)
        else
            x <- subsetByColData(x, j)
    }
    if (!missing(k)) {
        x <- subsetByAssay(x, k)
    }
    if (!missing(i)) {
        x <- subsetByRow(x, y = i, ...)
    }
    if (drop) {
        x <- .dropEmptyMASE(x, warn = TRUE)
    }
    x
}

#' @exportMethod [
setMethod("[", c("MultiAssaySpatialExperiment", "ANY", "ANY", "ANY"),
    .subsetMultiAssaySpatialExperiment)
