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
#' mapping (\code{assay}, \code{colname}, \code{element_type}, \code{region},
#' \code{instance_id}) linking experiment columns to rows in the spatial point
#' or shape layers via explicit foreign key.
#'
#' \strong{Terminology:} a \emph{specimen} is a row in \code{colData}
#' (\code{sampleMap$primary}); an \emph{observation} is a column in an assay
#' (\code{sampleMap$colname}).
#'
#' @section Slots:
#' \describe{
#'   \item{\code{ExperimentList}, \code{colData}, \code{sampleMap}, \code{drops}}{
#'     Inherited from \linkS4class{MultiAssayExperiment}.
#'   }
#'   \item{\code{images}}{
#'     A \linkS4class{RasterLayerList} of user-attached image rasters (not the
#'     same as \code{imgData}, which stores specimen-level image metadata from
#'     vendor readers).
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
#'     Optional \code{DataFrame} linking specimens to images (same structure as
#'     \code{SpatialExperiment}, including \code{sample_id}, \code{image_id}, and
#'     \code{scaleFactor}). Populated by vendor readers; distinct from
#'     \code{spatialImages()} raster layers.
#'   }
#'   \item{\code{spatialMap}}{
#'     Optional \code{DataFrame} for instance-level mapping with required columns:
#'     \describe{
#'       \item{\code{assay}}{Experiment name (factor; must be in \code{names(experiments)})}
#'       \item{\code{colname}}{Column identifier within that experiment}
#'       \item{\code{element_type}}{Spatial slot discriminator; one of \code{"points"} or \code{"shapes"}}
#'       \item{\code{region}}{Layer name within that element type (e.g., \code{"cells"}, \code{"nuclei"})}
#'       \item{\code{instance_id}}{Row identifier in that layer (character or integer; no NAs)}
#'     }
#'     The \code{element_type} column disambiguates spatial layers with the same
#'     name in different slots (e.g., \code{points$cells} vs. \code{shapes$cells}).
#'     This explicit routing eliminates ambiguity and enables robust foreign key
#'     validation.
#'   }
#' }
#'
#' @section Terminology and spatialdata alignment:
#' MASE adopts terminology from Python spatialdata where practical while following
#' Bioconductor S4 conventions:
#' \itemize{
#'   \item \code{element_type}: Matches spatialdata's element type discriminator
#'     (values: \code{"points"}, \code{"shapes"}, \code{"images"}, \code{"labels"}).
#'     Currently only \code{"points"} and \code{"shapes"} are supported in
#'     \code{spatialMap} linkage; images and labels are reference layers accessed
#'     via \code{spatialImages()} and \code{spatialLabels()}.
#'   \item \code{region}: Corresponds to spatialdata's region (layer name within an
#'     element type).
#'   \item \code{instance_id}: Matches spatialdata's instance identifier (row ID
#'     within a region). Recommended types: integer, character, or factor. Must be
#'     unique within each (\code{element_type}, \code{region}) combination.
#'   \item Class naming: MASE uses "Layer" suffix (\code{PointsLayerList},
#'     \code{ShapesLayerList}) rather than "Element" to avoid R naming collisions.
#'     See \code{?SpatialLayerList} for details.
#'   \item Accessor naming: MASE uses \code{spatialPoints()}, \code{spatialShapes()}
#'     (with "spatial" prefix) as S4 generics, while spatialdata uses bare attribute
#'     access (\code{sdata.points}, \code{sdata.shapes}). This follows Bioconductor
#'     conventions and avoids conflicts with base R functions.
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
#'         A \linkS4class{DataFrame} of specimen-level metadata (one row per
#'         primary identifier in \code{sampleMap$primary}).
#'       }
#'       \item{\code{sampleMap}}{
#'         A \linkS4class{DataFrame} with columns \code{assay}, \code{primary},
#'         \code{colname}.
#'       }
#'       \item{\code{metadata}, \code{drops}}{
#'         Additional \linkS4class{MultiAssayExperiment} arguments.
#'       }
#'       \item{\code{images}, \code{labels}}{
#'         \linkS4class{RasterLayerList} of user-attached raster layers (distinct
#'         from \code{imgData}, which vendor readers use for image metadata).
#'       }
#'       \item{\code{points}, \code{shapes}}{
#'         \linkS4class{PointsLayerList} and \linkS4class{ShapesLayerList}
#'         of spatial layers.
#'       }
#'       \item{\code{imgData}}{
#'         Optional \code{DataFrame} linking specimens to images (populated by
#'         vendor readers; use \code{getImg()} for SPE-compatible access).
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
#'     Get or set the optional \code{DataFrame} linking specimens to images
#'     (\code{SpatialExperiment}-compatible columns such as \code{sample_id},
#'     \code{image_id}, and \code{scaleFactor}).
#'   }
#'   \item{\code{spatialMap(x)}, \code{spatialMap(x) <- value}:}{
#'     Get or set the optional \code{DataFrame} for instance-level mapping.
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
#' mat <- matrix(rnorm(20), 5, 4,
#'     dimnames = list(paste0("G", 1:5), paste0("obs", 1:4)))
#' pts <- S4Vectors::DataFrame(
#'     x = 1:4, y = 1:4, instance_id = paste0("obs", 1:4))
#' mase <- MultiAssaySpatialExperiment(
#'     experiments = ExperimentList(rna = mat),
#'     colData = S4Vectors::DataFrame(
#'         tissue = c("core", "margin"),
#'         row.names = c("specimen_A", "specimen_B")),
#'     sampleMap = S4Vectors::DataFrame(
#'         assay = factor(rep("rna", 4), "rna"),
#'         primary = rep(c("specimen_A", "specimen_B"), each = 2),
#'         colname = paste0("obs", 1:4)),
#'     points = PointsLayerList(coords = pts),
#'     spatialMap = S4Vectors::DataFrame(
#'         assay = factor(rep("rna", 4), "rna"),
#'         colname = paste0("obs", 1:4),
#'         element_type = "points",
#'         region = "coords",
#'         instance_id = paste0("obs", 1:4))
#' )
#' mase
#' spatialPointNames(mase)
#' mase[, "specimen_A"]
#'
#' @return
#' The \code{MultiAssaySpatialExperiment()} constructor returns a
#' \linkS4class{MultiAssaySpatialExperiment}. The layer accessors return the
#' corresponding container: \code{spatialImages()} / \code{spatialLabels()} a
#' \linkS4class{RasterLayerList}, \code{spatialPoints()} a
#' \linkS4class{PointsLayerList}, and \code{spatialShapes()} a
#' \linkS4class{ShapesLayerList}; the single-layer accessors
#' (\code{spatialImage()}, \code{spatialPoint()}, and so on) return one layer,
#' and the \code{*Names()} accessors a character vector. \code{imgData()} and
#' \code{spatialMap()} return a \code{DataFrame} or \code{NULL}. Replacement
#' methods (the \code{<-} forms) return the updated object, and \code{show()}
#' returns \code{NULL} invisibly.
#'
#' @aliases MultiAssaySpatialExperiment-class
#'
#' @aliases MultiAssaySpatialExperiment
#'
#' @aliases spatialImages,MultiAssaySpatialExperiment-method
#' @aliases spatialImages<-,MultiAssaySpatialExperiment-method
#' @aliases spatialImage,MultiAssaySpatialExperiment-method
#' @aliases spatialImage<-,MultiAssaySpatialExperiment-method
#' @aliases spatialImageNames,MultiAssaySpatialExperiment-method
#'
#' @aliases spatialLabels,MultiAssaySpatialExperiment-method
#' @aliases spatialLabels<-,MultiAssaySpatialExperiment-method
#' @aliases spatialLabel,MultiAssaySpatialExperiment-method
#' @aliases spatialLabel<-,MultiAssaySpatialExperiment-method
#' @aliases spatialLabelNames,MultiAssaySpatialExperiment-method
#'
#' @aliases spatialPoints,MultiAssaySpatialExperiment-method
#' @aliases spatialPoints<-,MultiAssaySpatialExperiment-method
#' @aliases spatialPoint,MultiAssaySpatialExperiment-method
#' @aliases spatialPoint<-,MultiAssaySpatialExperiment-method
#' @aliases spatialPointNames,MultiAssaySpatialExperiment-method
#'
#' @aliases spatialShapes,MultiAssaySpatialExperiment-method
#' @aliases spatialShapes<-,MultiAssaySpatialExperiment-method
#' @aliases spatialShape,MultiAssaySpatialExperiment-method
#' @aliases spatialShape<-,MultiAssaySpatialExperiment-method
#' @aliases spatialShapeNames,MultiAssaySpatialExperiment-method
#'
#' @aliases imgData,MultiAssaySpatialExperiment-method
#' @aliases imgData<-,MultiAssaySpatialExperiment,ANY-method
#'
#' @aliases spatialMap,MultiAssaySpatialExperiment-method
#' @aliases spatialMap<-,MultiAssaySpatialExperiment-method
#'
#' @aliases show,MultiAssaySpatialExperiment-method
#'
#' @seealso
#' \itemize{
#'   \item \link{MultiAssaySpatialExperiment-subset} for subsetting methods
#'   \item \link{MultiAssaySpatialExperiment-combine} for combining methods
#'   \item \linkS4class{RasterLayerList}, \linkS4class{PointsLayerList},
#'         \linkS4class{ShapesLayerList} for spatial layer containers
#' }
#'
#' @keywords classes methods
#'
#' @name MultiAssaySpatialExperiment-class
NULL

#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @importClassesFrom S4Vectors DataFrame_OR_NULL
#' @importFrom S4Vectors SimpleList
#' @exportClass MultiAssaySpatialExperiment
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

## SPATIALMAP
## 1.i. spatialMap must have required columns
.checkSpatialMapColumns <- function(object) {
    errors <- character()
    spmap <- spatialMap(object)

    if (!is.null(spmap) && nrow(spmap) > 0L) {
        if (!all(SPATIAL_MAP_COLS %in% colnames(spmap))) {
            msg <- sprintf(
                "'spatialMap' must have columns: %s",
                paste(SPATIAL_MAP_COLS, collapse = ", ")
            )
            errors <- c(errors, msg)
        }
    }

    if (!length(errors)) NULL else errors
}

## 1.ii. spatialMap (assay, colname) must be found in sampleMap
.checkSpatialMapFK <- function(object) {
    errors <- character()
    spmap <- spatialMap(object)

    if (!is.null(spmap) && nrow(spmap) > 0L &&
        all(c("assay", "colname") %in% colnames(spmap))) {
        sm <- sampleMap(object)
        if (nrow(sm) > 0L) {
            sm_key <- paste(sm[["assay"]], sm[["colname"]], sep = "\r")
            linked <- !is.na(spmap[["colname"]])
            sp_key <- paste(spmap[["assay"]], spmap[["colname"]], sep = "\r")
            bad <- linked & !sp_key %in% sm_key
            if (any(bad)) {
                msg <- paste(
                    "spatialMap (assay, colname) must be in sampleMap;",
                    "found", sum(bad), "row(s) not in sampleMap"
                )
                errors <- c(errors, msg)
            }
        }
    }

    if (!length(errors)) NULL else errors
}

## 1.iii. spatialMap element_type must be valid
.checkElementTypes <- function(object) {
    errors <- character()
    spmap <- spatialMap(object)

    if (!is.null(spmap) && nrow(spmap) > 0L) {
        if (!"element_type" %in% colnames(spmap)) {
            errors <- c(errors, "'spatialMap' must have 'element_type' column")
        } else {
            invalid_types <- !spmap[["element_type"]] %in% VALID_ELEMENT_TYPES
            if (any(invalid_types)) {
                msg <- sprintf(
                    "spatialMap$element_type contains invalid value(s): %s; must be one of: %s",
                    paste(unique(spmap[["element_type"]][invalid_types]), collapse = ", "),
                    paste(VALID_ELEMENT_TYPES, collapse = ", ")
                )
                errors <- c(errors, msg)
            }
        }
    }

    if (!length(errors)) NULL else errors
}

## 1.iv. spatialMap instance_id must not contain NAs and must reference valid
## spatial layer rows
.checkInstanceIds <- function(object) {
    errors <- character()
    spmap <- spatialMap(object)

    if (!is.null(spmap) && nrow(spmap) > 0L &&
        all(c("element_type", "region", "instance_id") %in% colnames(spmap))) {

        ## Check NAs in spatialMap instance_id
        if (anyNA(spmap[["instance_id"]])) {
            errors <- c(errors, "spatialMap$instance_id contains NA values")
        }

        ## Validate (element_type, region, instance_id) foreign keys
        ## Route by element_type to avoid ambiguity
        for (et in unique(spmap[["element_type"]])) {
            if (!et %in% VALID_ELEMENT_TYPES) next

            et_rows <- spmap[["element_type"]] == et
            et_spmap <- spmap[et_rows, , drop = FALSE]

            ## Get the spatial slot for this element type
            slot_data <- switch(et,
                points = spatialPoints(object),
                shapes = spatialShapes(object),
                NULL  ## Extensibility: add "labels" = spatialLabels(object) here
            )

            if (is.null(slot_data)) next
            slot_names <- names(slot_data)

            ## Check each region in this element type
            for (reg in unique(et_spmap[["region"]])) {
                ## Region must exist in slot
                if (!reg %in% slot_names) {
                    msg <- sprintf(
                        "spatialMap: region '%s' (element_type = '%s') not found in %s slot",
                        reg, et, et
                    )
                    errors <- c(errors, msg)
                    next
                }

                ## Get layer and validate instance_id values
                layer <- slot_data[[reg]]
                if (is.null(layer) || nrow(layer) == 0L) next

                ## Check for NAs in layer instance_id column (if present)
                if ("instance_id" %in% colnames(layer)) {
                    if (anyNA(layer[["instance_id"]])) {
                        msg <- sprintf(
                            "instance_id column in %s$%s contains NA values",
                            et, reg
                        )
                        errors <- c(errors, msg)
                    }

                    ## Warn if instance_id type is unusual (informational only)
                    id_col <- layer[["instance_id"]]
                    id_type <- class(id_col)[1L]
                    recommended_types <- c("integer", "numeric", "character", 
                                          "factor", "Rle")
                    if (!id_type %in% recommended_types) {
                        warning(wmsg(
                            "instance_id in ", et, "$", reg, 
                            " has unusual type: ", id_type, "; ",
                            "recommend integer, character, or factor for ",
                            "spatialdata interoperability"
                        ), call. = FALSE)
                    }
                }

                ## Check instance_id values exist in layer
                region_rows <- et_spmap[["region"]] == reg
                region_instances <- unique(et_spmap[["instance_id"]][region_rows])

                ## Get layer instance IDs (from instance_id column or row numbers)
                layer_ids <- if ("instance_id" %in% colnames(layer)) {
                    as.character(layer[["instance_id"]])
                } else {
                    as.character(seq_len(nrow(layer)))
                }

                ## Check all spatialMap instance_id values exist in layer
                missing <- !as.character(region_instances) %in% layer_ids
                if (any(missing)) {
                    msg <- sprintf(
                        "spatialMap: instance_id value(s) not found in %s$%s: %s",
                        et, reg,
                        paste(as.character(region_instances[missing]), collapse = ", ")
                    )
                    errors <- c(errors, msg)
                }
            }
        }
    }

    if (!length(errors)) NULL else errors
}

## IMGDATA
## 2.i. imgData sample_id must be found in colData rownames
.checkImgData <- function(object) {
    errors <- character()
    img <- imgData(object)

    if (!is.null(img) && nrow(img) > 0L && "sample_id" %in% colnames(img)) {
        cd <- colData(object)
        primaries <- if ("sample_id" %in% colnames(cd)) {
            unique(as.character(cd[["sample_id"]]))
        } else {
            rownames(cd)
        }
        bad <- !img[["sample_id"]] %in% primaries
        if (any(bad) && length(primaries) > 0L) {
            msg <- paste(
                "imgData sample_id must be in colData rownames (primary_id);",
                "found", sum(bad), "row(s) not in colData"
            )
            errors <- c(errors, msg)
        }
    }

    if (!length(errors)) NULL else errors
}

#' @importFrom MultiAssayExperiment sampleMap
#' @importFrom SpatialExperiment imgData
#' @importFrom SummarizedExperiment colData
.validMultiAssaySpatialExperiment <- function(object) {
    c(.checkSpatialMapColumns(object),
      .checkSpatialMapFK(object),
      .checkElementTypes(object),
      .checkInstanceIds(object),
      .checkImgData(object))
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
    x <- replaceSlots(x, spatialMap = value, check = FALSE)
    validObject(x)
    x
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
