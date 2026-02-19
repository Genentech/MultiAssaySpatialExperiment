### =========================================================================
### MultiAssaySpatialExperiment class
### -------------------------------------------------------------------------
###
### Extends MultiAssayExperiment with spatial elements: images, labels,
### points, shapes, imgData, spatialMap. Central mapping tables for RI.
###
### -------------------------------------------------------------------------

#' @include SpatialGenerics.R
#' @include SpatialLayerList-class.R
NULL

replaceSlots <- BiocGenerics:::replaceSlots

#' MultiAssaySpatialExperiment objects
#'
#' @description
#' The MultiAssaySpatialExperiment class extends \linkS4class{MultiAssayExperiment}
#' with spatial context: raster data (images, labels), points, and shapes.
#'
#' @details
#' MultiAssaySpatialExperiment provides slots for spatial layers that complement
#' the \code{ExperimentList}. The \code{sampleMap} links experiment columns to
#' primary identifiers; the optional \code{spatialMap} provides instance-level
#' mapping (\code{assay}, \code{colname}, \code{region}, \code{instance_id})
#' linking experiment columns to rows in the spatial point or shape layers.
#'
#' @section Slots:
#' \describe{
#'   \item{\code{ExperimentList}, \code{colData}, \code{sampleMap}, \code{drops}}{
#'     Inherited from \linkS4class{MultiAssayExperiment}.
#'   }
#'   \item{\code{images}}{
#'     A \linkS4class{RasterLayerList} of image layers.
#'   }
#'   \item{\code{labels}}{
#'     A \linkS4class{RasterLayerList} of label/mask layers.
#'   }
#'   \item{\code{points}}{
#'     A \linkS4class{PointsLayerList} of point layers.
#'   }
#'   \item{\code{shapes}}{
#'     A \linkS4class{ShapesLayerList} of shape layers.
#'   }
#'   \item{\code{imgData}}{
#'     Optional \code{DataFrame} of sample-linked images (same structure as
#'     \code{SpatialExperiment}); for compatibility when experiments are
#'     SpatialExperiments.
#'   }
#'   \item{\code{spatialMap}}{
#'     Optional \code{DataFrame} for instance-level mapping:
#'     (\code{assay}, \code{colname}, \code{region}, \code{instance_id})
#'     linking experiment columns to spatial layer rows.
#'   }
#' }
#'
#' @section Constructor:
#' \describe{
#'   \item{\code{MultiAssaySpatialExperiment(experiments = ExperimentList(),
#'         colData = DataFrame(), sampleMap = DataFrame(...), metadata = list(),
#'         drops = list(), images = RasterLayerList(), labels = RasterLayerList(),
#'         points = PointsLayerList(), shapes = ShapesLayerList(),
#'         imgData = NULL, spatialMap = NULL)}:}{
#'     Creates a MultiAssaySpatialExperiment object.
#'     \describe{
#'       \item{\code{experiments}}{
#'         An \linkS4class{ExperimentList} or list of experiments.
#'       }
#'       \item{\code{colData}}{
#'         A \linkS4class{DataFrame} of subject-level metadata.
#'       }
#'       \item{\code{sampleMap}}{
#'         A \linkS4class{DataFrame} with columns \code{assay}, \code{primary},
#'         \code{colname}.
#'       }
#'       \item{\code{metadata}, \code{drops}}{
#'         Additional \linkS4class{MultiAssayExperiment} arguments.
#'       }
#'       \item{\code{images}, \code{labels}}{
#'         \linkS4class{RasterLayerList} of image and label layers.
#'       }
#'       \item{\code{points}, \code{shapes}}{
#'         \linkS4class{PointsLayerList} and \linkS4class{ShapesLayerList}
#'         of spatial layers.
#'       }
#'       \item{\code{imgData}}{
#'         Optional \code{DataFrame} for sample-linked images.
#'       }
#'       \item{\code{spatialMap}}{
#'         Optional \code{DataFrame} for instance-level mapping.
#'       }
#'     }
#'   }
#' }
#'
#' @section Accessors:
#' In the code snippets below, \code{x} is a MultiAssaySpatialExperiment object:
#' \describe{
#'   \item{\code{spatialImages(x)}, \code{spatialImages(x) <- value}:}{
#'     Get or set the full \linkS4class{RasterLayerList} of image layers.
#'   }
#'   \item{\code{spatialLabels(x)}, \code{spatialLabels(x) <- value}:}{
#'     Get or set the full \linkS4class{RasterLayerList} of label layers.
#'   }
#'   \item{\code{spatialPoints(x)}, \code{spatialPoints(x) <- value}:}{
#'     Get or set the full \linkS4class{PointsLayerList} of point layers.
#'   }
#'   \item{\code{spatialShapes(x)}, \code{spatialShapes(x) <- value}:}{
#'     Get or set the full \linkS4class{ShapesLayerList} of shape layers.
#'   }
#'   \item{\code{spatialImage(x, i)}, \code{spatialImage(x, i) <- value}:}{
#'     Get or set a single image layer by index \code{i} (character or integer).
#'   }
#'   \item{\code{spatialLabel(x, i)}, \code{spatialLabel(x, i) <- value}:}{
#'     Get or set a single label layer by index \code{i}.
#'   }
#'   \item{\code{spatialPoint(x, i)}, \code{spatialPoint(x, i) <- value}:}{
#'     Get or set a single point layer by index \code{i}.
#'   }
#'   \item{\code{spatialShape(x, i)}, \code{spatialShape(x, i) <- value}:}{
#'     Get or set a single shape layer by index \code{i}.
#'   }
#'   \item{\code{spatialImageNames(x)}, \code{spatialLabelNames(x)},
#'     \code{spatialPointNames(x)}, \code{spatialShapeNames(x)}:}{
#'     Get the names of the image, label, point, and shape layers, respectively.
#'   }
#'   \item{\code{imgData(x)}, \code{imgData(x) <- value}:}{
#'     Get or set the optional \code{DataFrame} of sample-linked images.
#'   }
#'   \item{\code{spatialMap(x)}, \code{spatialMap(x) <- value}:}{
#'     Get or set the optional \code{DataFrame} for instance-level mapping.
#'   }
#' }
#'
#' @section Combining:
#' \describe{
#'   \item{\code{c(x, ..., sampleMap = NULL, mapFrom = NULL)}:}{
#'     When \code{...} is a single MultiAssaySpatialExperiment, merges the two
#'     objects (experiments, colData, sampleMap, spatial layers, imgData,
#'     spatialMap). Experiment names must be unique across objects.
#'     Otherwise delegates to \code{\link[MultiAssayExperiment]{c,MultiAssayExperiment-method}}.
#'   }
#' }
#'
#' @section Displaying:
#' The \code{show()} method prints the inherited MultiAssayExperiment summary
#' followed by counts of spatial images, labels, points, shapes, and whether
#' \code{imgData} and \code{spatialMap} are present.
#'
#' @author Patrick Aboyoun
#'
#' @examples
#' mat <- matrix(rnorm(20), 5, 4, dimnames = list(paste0("G", 1:5), paste0("S", 1:4)))
#' pts <- S4Vectors::DataFrame(x = c(1,2,3,4), y = c(1,2,3,4), instance_id = paste0("S", 1:4))
#' mase <- MultiAssaySpatialExperiment(
#'     experiments = ExperimentList(assay1 = mat),
#'     colData = S4Vectors::DataFrame(row.names = paste0("S", 1:4)),
#'     sampleMap = S4Vectors::DataFrame(
#'         assay = factor("assay1", "assay1"),
#'         primary = paste0("S", 1:4),
#'         colname = paste0("S", 1:4)
#'     ),
#'     points = PointsLayerList(coords = pts)
#' )
#' mase
#' spatialPointNames(mase)
#' mase[, c("S1", "S2")]
#'
#' @aliases
#' MultiAssaySpatialExperiment-class
#'
#' MultiAssaySpatialExperiment
#'
#' spatialImages,MultiAssaySpatialExperiment-method
#' spatialImages<-,MultiAssaySpatialExperiment-method
#' spatialImage,MultiAssaySpatialExperiment-method
#' spatialImage<-,MultiAssaySpatialExperiment-method
#' spatialImageNames,MultiAssaySpatialExperiment-method
#'
#' spatialLabels,MultiAssaySpatialExperiment-method
#' spatialLabels<-,MultiAssaySpatialExperiment-method
#' spatialLabel,MultiAssaySpatialExperiment-method
#' spatialLabel<-,MultiAssaySpatialExperiment-method
#' spatialLabelNames,MultiAssaySpatialExperiment-method
#'
#' spatialPoints,MultiAssaySpatialExperiment-method
#' spatialPoints<-,MultiAssaySpatialExperiment-method
#' spatialPoint,MultiAssaySpatialExperiment-method
#' spatialPoint<-,MultiAssaySpatialExperiment-method
#' spatialPointNames,MultiAssaySpatialExperiment-method
#'
#' spatialShapes,MultiAssaySpatialExperiment-method
#' spatialShapes<-,MultiAssaySpatialExperiment-method
#' spatialShape,MultiAssaySpatialExperiment-method
#' spatialShape<-,MultiAssaySpatialExperiment-method
#' spatialShapeNames,MultiAssaySpatialExperiment-method
#'
#' imgData,MultiAssaySpatialExperiment-method
#' imgData<-,MultiAssaySpatialExperiment,ANY-method
#'
#' spatialMap,MultiAssaySpatialExperiment-method
#' spatialMap<-,MultiAssaySpatialExperiment-method
#'
#' c,MultiAssaySpatialExperiment-method
#'
#' show,MultiAssaySpatialExperiment-method
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[MultiAssayExperiment]{MultiAssayExperiment}} for the base class
#'   \item \link{MultiAssaySpatialExperiment-subset} for subsetting methods
#'   \item \linkS4class{RasterLayerList}, \linkS4class{PointsLayerList},
#'         \linkS4class{ShapesLayerList} for spatial layer containers
#' }
#'
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @importClassesFrom S4Vectors DataFrame_OR_NULL
#'
#' @keywords classes methods
#'
#' @exportClass MultiAssaySpatialExperiment
#'
#' @name MultiAssaySpatialExperiment-class
NULL

setClass("MultiAssaySpatialExperiment", contains = "MultiAssayExperiment",
    slots = c(images   = "RasterLayerList",
              labels   = "RasterLayerList",
              points   = "PointsLayerList",
              shapes   = "ShapesLayerList",
              imgData  = "DataFrame_OR_NULL",
              spatialMap = "DataFrame_OR_NULL"),
    prototype = prototype(
        images = new("RasterLayerList",
                     SimpleList(structure(list(), names = character()))),
        labels = new("RasterLayerList",
                     SimpleList(structure(list(), names = character()))),
        points = new("PointsLayerList",
                     SimpleList(structure(list(), names = character()))),
        shapes = new("ShapesLayerList",
                     SimpleList(structure(list(), names = character()))),
        imgData = NULL,
        spatialMap = NULL
    )
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

SPATIAL_MAP_COLS <- c("assay", "colname", "region", "instance_id")

#' @importFrom MultiAssayExperiment sampleMap
#' @importFrom SpatialExperiment imgData
#' @importFrom SummarizedExperiment colData
.validMultiAssaySpatialExperiment <- function(object) {
    msg <- NULL

    spmap <- spatialMap(object)
    if (!is.null(spmap) && nrow(spmap) > 0L) {
        if (!all(SPATIAL_MAP_COLS %in% colnames(spmap))) {
            msg <- c(msg, sprintf(
                "'spatialMap' must have columns: %s",
                paste(SPATIAL_MAP_COLS, collapse = ", ")
            ))
        } else {
            sm <- sampleMap(object)
            if (nrow(sm) > 0L) {
                sm_key <- paste(sm[["assay"]], sm[["colname"]], sep = "\r")
                sp_key <- paste(spmap[["assay"]], spmap[["colname"]], sep = "\r")
                bad <- !sp_key %in% sm_key
                if (any(bad)) {
                    msg <- c(msg, paste(
                        "spatialMap (assay, colname) must be in sampleMap;",
                        "found", sum(bad), "row(s) not in sampleMap"
                    ))
                }
            }
            pts <- spatialPoints(object)
            shps <- spatialShapes(object)
            pt_names <- names(pts)
            shp_names <- names(shps)
            valid_regions <- c(pt_names, shp_names)
            for (i in seq_len(nrow(spmap))) {
                reg <- spmap[["region"]][[i]]
                if (!reg %in% valid_regions) next
                inst <- spmap[["instance_id"]][[i]]
                if (reg %in% pt_names) {
                    el <- pts[[reg]]
                    if (!is.null(el) && nrow(el) > 0L) {
                        id_col <- if ("instance_id" %in% colnames(el))
                            el[["instance_id"]] else seq_len(nrow(el))
                        id_chr <- as.character(id_col)
                        if (!as.character(inst) %in% id_chr)
                            msg <- c(msg, sprintf(
                                "spatialMap row %d: instance_id %s not in points[[%s]]",
                                i, as.character(inst), reg
                            ))
                    }
                } else {
                    el <- shps[[reg]]
                    if (!is.null(el) && nrow(el) > 0L) {
                        id_col <- if ("instance_id" %in% colnames(el))
                            el[["instance_id"]] else seq_len(nrow(el))
                        id_chr <- as.character(id_col)
                        if (!as.character(inst) %in% id_chr)
                            msg <- c(msg, sprintf(
                                "spatialMap row %d: instance_id %s not in shapes[[%s]]",
                                i, as.character(inst), reg
                            ))
                    }
                }
            }
        }
    }

    img <- imgData(object)
    if (!is.null(img) && nrow(img) > 0L && "sample_id" %in% colnames(img)) {
        primaries <- rownames(colData(object))
        bad <- !img[["sample_id"]] %in% primaries
        if (any(bad) && length(primaries) > 0L) {
            msg <- c(msg, paste(
                "imgData sample_id must be in colData rownames (primary_id);",
                "found", sum(bad), "row(s) not in colData"
            ))
        }
    }

    msg
}

#' @importFrom S4Vectors setValidity2
setValidity2("MultiAssaySpatialExperiment", .validMultiAssaySpatialExperiment)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

#' @importFrom MultiAssayExperiment ExperimentList MultiAssayExperiment
#' @importFrom S4Vectors DataFrame
#' @export
MultiAssaySpatialExperiment <-
function(experiments = ExperimentList(),
         colData = DataFrame(),
         sampleMap = DataFrame(assay = factor(),
                               primary = character(),
                               colname = character()),
         metadata = list(),
         drops = list(),
         images = RasterLayerList(),
         labels = RasterLayerList(),
         points = PointsLayerList(),
         shapes = ShapesLayerList(),
         imgData = NULL,
         spatialMap = NULL)
{
    mae <- MultiAssayExperiment(experiments = experiments,
                                colData = colData,
                                sampleMap = sampleMap,
                                metadata = metadata,
                                drops = drops)
    new("MultiAssaySpatialExperiment",
        mae,
        images = as(images, "RasterLayerList"),
        labels = as(labels, "RasterLayerList"),
        points = as(points, "PointsLayerList"),
        shapes = as(shapes, "ShapesLayerList"),
        imgData = imgData,
        spatialMap = spatialMap
    )
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

### spatialImages / spatialImage

#' @export
setMethod("spatialImages", "MultiAssaySpatialExperiment", function(x) slot(x, "images"))

#' @export
setReplaceMethod("spatialImages", "MultiAssaySpatialExperiment", function(x, value) {
    replaceSlots(x, images = as(value, "RasterLayerList"), check = FALSE)
})

#' @export
setMethod("spatialImage", c("MultiAssaySpatialExperiment", "ANY"), function(x, i) spatialImages(x)[[i]])

#' @export
setReplaceMethod("spatialImage", c("MultiAssaySpatialExperiment", "ANY"), function(x, i, value) {
    imgs <- spatialImages(x)
    imgs[[i]] <- value
    spatialImages(x) <- imgs
    x
})

#' @export
setMethod("spatialImageNames", "MultiAssaySpatialExperiment", function(x) names(spatialImages(x)))

### spatialLabels / spatialLabel

#' @export
setMethod("spatialLabels", "MultiAssaySpatialExperiment", function(x) slot(x, "labels"))

#' @export
setReplaceMethod("spatialLabels", "MultiAssaySpatialExperiment", function(x, value) {
    replaceSlots(x, labels = as(value, "RasterLayerList"), check = FALSE)
})

#' @export
setMethod("spatialLabel", c("MultiAssaySpatialExperiment", "ANY"), function(x, i) spatialLabels(x)[[i]])

#' @export
setMethod("spatialLabel", c("MultiAssaySpatialExperiment", "ANY"), function(x, i) spatialLabels(x)[[i]])

#' @export
setReplaceMethod("spatialLabel", c("MultiAssaySpatialExperiment", "ANY"), function(x, i, value) {
    lbls <- spatialLabels(x)
    lbls[[i]] <- value
    spatialLabels(x) <- lbls
    x
})

#' @export
setMethod("spatialLabelNames", "MultiAssaySpatialExperiment", function(x) names(spatialLabels(x)))

### spatialPoints / spatialPoint

#' @export
setMethod("spatialPoints", "MultiAssaySpatialExperiment", function(x) slot(x, "points"))

#' @export
setReplaceMethod("spatialPoints", "MultiAssaySpatialExperiment", function(x, value) {
    replaceSlots(x, points = as(value, "PointsLayerList"), check = FALSE)
})

#' @export
setMethod("spatialPoint", c("MultiAssaySpatialExperiment", "ANY"), function(x, i) spatialPoints(x)[[i]])

#' @export
setReplaceMethod("spatialPoint", c("MultiAssaySpatialExperiment", "ANY"), function(x, i, value) {
    pts <- spatialPoints(x)
    pts[[i]] <- value
    spatialPoints(x) <- pts
    x
})

#' @export
setMethod("spatialPointNames", "MultiAssaySpatialExperiment", function(x) names(spatialPoints(x)))

### spatialShapes / spatialShape

#' @export
setMethod("spatialShapes", "MultiAssaySpatialExperiment", function(x) slot(x, "shapes"))

#' @export
setReplaceMethod("spatialShapes", "MultiAssaySpatialExperiment", function(x, value) {
    replaceSlots(x, shapes = as(value, "ShapesLayerList"), check = FALSE)
})

#' @export
setMethod("spatialShape", c("MultiAssaySpatialExperiment", "ANY"), function(x, i) spatialShapes(x)[[i]])

#' @export
setReplaceMethod("spatialShape", c("MultiAssaySpatialExperiment", "ANY"), function(x, i, value) {
    shps <- spatialShapes(x)
    shps[[i]] <- value
    spatialShapes(x) <- shps
    x
})

#' @export
setMethod("spatialShapeNames", "MultiAssaySpatialExperiment", function(x) names(spatialShapes(x)))

### imgData

#' @importFrom SpatialExperiment imgData
#' @export
setMethod("imgData", "MultiAssaySpatialExperiment", function(x) slot(x, "imgData"))

#' @export
setReplaceMethod("imgData", "MultiAssaySpatialExperiment", function(x, value) {
    replaceSlots(x, imgData = value, check = FALSE)
})

### spatialMap

#' @export
setMethod("spatialMap", "MultiAssaySpatialExperiment", function(x) slot(x, "spatialMap"))

#' @export
setReplaceMethod("spatialMap", "MultiAssaySpatialExperiment", function(x, value) {
    replaceSlots(x, spatialMap = value, check = FALSE)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### c concatenation
###

#' @importFrom MultiAssayExperiment experiments sampleMap
#' @importFrom S4Vectors DataFrame metadata
#' @importFrom SpatialExperiment imgData
#' @importFrom SummarizedExperiment colData
.mergeMASE <- function(x, y) {
    if (!is(y, "MultiAssaySpatialExperiment"))
        stop("'c' with two MultiAssaySpatialExperiment requires both to be ",
             "MultiAssaySpatialExperiment")
    if (any(names(x) %in% names(y)))
        stop("Provide unique experiment names")
    expz <- c(experiments(x), experiments(y))
    sampz <- rbind(sampleMap(x), sampleMap(y))
    coldx <- colData(x)
    coldy <- colData(y)
    cdatz <- merge(x = coldx, y = coldy,
        by = c("row.names", intersect(names(coldx), names(coldy))),
        all = TRUE, sort = FALSE, stringsAsFactors = FALSE)
    rownames(cdatz) <- cdatz[["Row.names"]]
    cdatz <- DataFrame(cdatz[, names(cdatz) != "Row.names", drop = FALSE])
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
    new("MultiAssaySpatialExperiment",
        ExperimentList = expz,
        colData = cdatz,
        sampleMap = sampz,
        metadata = metaz,
        images = c(spatialImages(x), spatialImages(y)),
        labels = c(spatialLabels(x), spatialLabels(y)),
        points = c(spatialPoints(x), spatialPoints(y)),
        shapes = c(spatialShapes(x), spatialShapes(y)),
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
### Display
###

#' @export
setMethod("show", "MultiAssaySpatialExperiment", function(object) {
    callNextMethod()
    nimg <- length(spatialImages(object))
    nlbl <- length(spatialLabels(object))
    npt <- length(spatialPoints(object))
    nsh <- length(spatialShapes(object))
    cat("Spatial elements:\n",
        "  spatialImages: ", nimg, if (nimg != 1L) " elements\n" else " element\n",
        "  spatialLabels: ", nlbl, if (nlbl != 1L) " elements\n" else " element\n",
        "  spatialPoints: ", npt, if (npt != 1L) " elements\n" else " element\n",
        "  spatialShapes: ", nsh, if (nsh != 1L) " elements\n" else " element\n",
        "  imgData: ", if (is.null(imgData(object))) "NULL\n" else "present\n",
        "  spatialMap: ", if (is.null(spatialMap(object))) "NULL\n" else "present\n",
        sep = "")
})
