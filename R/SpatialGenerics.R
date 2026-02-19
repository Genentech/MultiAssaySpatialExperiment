### =========================================================================
### Generics for spatial accessors
### -------------------------------------------------------------------------
###
### spatialImages, spatialLabels, spatialPoints, spatialShapes, spatialMap,
### imgData (from SpatialExperiment), and single-element / name accessors.
###
### -------------------------------------------------------------------------

#' Access spatial element slots
#'
#' @description
#' Generic functions to access spatialImages, spatialLabels, spatialPoints,
#' spatialShapes, and spatialMap slots of
#' \linkS4class{MultiAssaySpatialExperiment}. The \code{imgData} generic
#' is imported from \pkg{SpatialExperiment}.
#'
#' @param x A \linkS4class{MultiAssaySpatialExperiment} object
#' @return For accessors: the slot or element contents. For replacement
#'   functions: the modified object.
#'
#' @details
#' Slot accessors (\code{spatialImages}, \code{spatialLabels}, etc.) return
#' the full element list. Single-element accessors (\code{spatialImage},
#' \code{spatialLabel}, etc.) take an index \code{i} (character or integer)
#' and return the element. Name accessors return character vectors of
#' element names. Replacement functions for slots take a full list; single-
#' element replacements take an index and the new value.
#'
#' For single-element accessors, a character index with no matching name
#' (e.g., \code{spatialImage(x, "nonexistent")}) returns \code{NULL}.
#' Integer indices that are out of bounds raise an error.
#'
#' @author Patrick Aboyoun
#'
#' @seealso
#' \linkS4class{MultiAssaySpatialExperiment},
#' \linkS4class{RasterLayerList},
#' \linkS4class{PointsLayerList},
#' \linkS4class{ShapesLayerList}
#'
#' @examples
#' mat <- matrix(rnorm(20), 5, 4, dimnames = list(paste0("G", 1:5), paste0("S", 1:4)))
#' expList <- ExperimentList(assay1 = mat)
#' cd <- DataFrame(row.names = paste0("S", 1:4))
#' sm <- DataFrame(assay = factor("assay1", "assay1"), primary = paste0("S", 1:4),
#'     colname = paste0("S", 1:4))
#' mase <- MultiAssaySpatialExperiment(
#'     experiments = expList,
#'     colData = cd,
#'     sampleMap = sm,
#'     images = RasterLayerList(brightfield = matrix(1:12, 3, 4))
#' )
#' spatialImageNames(mase)
#' spatialImage(mase, 1L)
#' spatialImage(mase, "brightfield")
#' spatialImages(mase) <- RasterLayerList(fluorescent = matrix(20:31, 3, 4))
#' spatialImageNames(mase)
#'
#' @import methods BiocGenerics
#' @importFrom SpatialExperiment imgData
#'
#' @name spatial-accessors
NULL

#' @rdname spatial-accessors
#' @export
setGeneric("spatialImages", function(x) standardGeneric("spatialImages"))

#' @rdname spatial-accessors
#' @export
setGeneric("spatialLabels", function(x) standardGeneric("spatialLabels"))

#' @rdname spatial-accessors
#' @export
setGeneric("spatialPoints", function(x) standardGeneric("spatialPoints"))

#' @rdname spatial-accessors
#' @export
setGeneric("spatialShapes", function(x) standardGeneric("spatialShapes"))

#' @rdname spatial-accessors
#' @export
setGeneric("spatialMap", function(x) standardGeneric("spatialMap"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Slot setters
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @rdname spatial-accessors
#' @param value Replacement value for the slot
#' @export
setGeneric("spatialImages<-", function(x, value) standardGeneric("spatialImages<-"))

#' @rdname spatial-accessors
#' @export
setGeneric("spatialLabels<-", function(x, value) standardGeneric("spatialLabels<-"))

#' @rdname spatial-accessors
#' @export
setGeneric("spatialPoints<-", function(x, value) standardGeneric("spatialPoints<-"))

#' @rdname spatial-accessors
#' @export
setGeneric("spatialShapes<-", function(x, value) standardGeneric("spatialShapes<-"))

#' @rdname spatial-accessors
#' @export
setGeneric("spatialMap<-", function(x, value) standardGeneric("spatialMap<-"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Single-element accessors and setters
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @rdname spatial-accessors
#' @param i Index (character or integer) for the element name
#' @export
setGeneric("spatialImage", function(x, i) standardGeneric("spatialImage"))

#' @rdname spatial-accessors
#' @export
setGeneric("spatialLabel", function(x, i) standardGeneric("spatialLabel"))

#' @rdname spatial-accessors
#' @export
setGeneric("spatialPoint", function(x, i) standardGeneric("spatialPoint"))

#' @rdname spatial-accessors
#' @export
setGeneric("spatialShape", function(x, i) standardGeneric("spatialShape"))

#' @rdname spatial-accessors
#' @export
setGeneric("spatialImageNames", function(x) standardGeneric("spatialImageNames"))

#' @rdname spatial-accessors
#' @export
setGeneric("spatialLabelNames", function(x) standardGeneric("spatialLabelNames"))

#' @rdname spatial-accessors
#' @export
setGeneric("spatialPointNames", function(x) standardGeneric("spatialPointNames"))

#' @rdname spatial-accessors
#' @export
setGeneric("spatialShapeNames", function(x) standardGeneric("spatialShapeNames"))

#' @rdname spatial-accessors
#' @param value Replacement value for the element
#' @export
setGeneric("spatialImage<-", function(x, i, value) standardGeneric("spatialImage<-"))

#' @rdname spatial-accessors
#' @export
setGeneric("spatialLabel<-", function(x, i, value) standardGeneric("spatialLabel<-"))

#' @rdname spatial-accessors
#' @export
setGeneric("spatialPoint<-", function(x, i, value) standardGeneric("spatialPoint<-"))

#' @rdname spatial-accessors
#' @export
setGeneric("spatialShape<-", function(x, i, value) standardGeneric("spatialShape<-"))
