### =========================================================================
### SpatialLayerList - Typed spatial element containers
### -------------------------------------------------------------------------
###
### RasterLayerList, PointsLayerList, ShapesLayerList for images,
### labels, points, and shapes. Class definitions only.
###
### -------------------------------------------------------------------------

#' @include SpatialGenerics.R
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Class definitions
###

#' Spatial layer list objects
#'
#' @description
#' The \linkS4class{SpatialLayerList} class is a virtual base extending
#' \linkS4class{SimpleList}. Concrete subclasses \linkS4class{RasterLayerList},
#' \linkS4class{PointsLayerList}, and \linkS4class{ShapesLayerList} provide
#' typed containers for raster data (images, labels), point coordinates, and
#' shape geometries.
#'
#' @details
#' All layers must be named and have unique names. Element-type validity is
#' enforced by each concrete class: \linkS4class{RasterLayerList} requires
#' \code{dim()} with length >= 2 (e.g., \code{matrix}, \code{array},
#' \code{DelayedArray}, \code{StoredSpatialImage}); \linkS4class{PointsLayerList}
#' requires \code{DataFrame} elements with at least 2 columns (typically x, y)
#' and often an instance identifier; \linkS4class{ShapesLayerList} requires
#' \code{DataFrame} elements (with \pkg{sf} a geometry column may be present).
#'
#' @section Constructors:
#' \describe{
#'   \item{\code{RasterLayerList(...)}:}{
#'     Creates a \linkS4class{RasterLayerList}. Accepts named arguments, a
#'     single named list, or a single \linkS4class{List}. Each element must
#'     have \code{dim()} with length >= 2.
#'   }
#'   \item{\code{PointsLayerList(...)}:}{
#'     Creates a \linkS4class{PointsLayerList}. Accepts named \code{DataFrame}
#'     elements or a single list/List. Each element must have at least 2 columns.
#'   }
#'   \item{\code{ShapesLayerList(...)}:}{
#'     Creates a \linkS4class{ShapesLayerList}. Accepts named \code{DataFrame}
#'     or \code{data.frame} elements (coerced to \code{DataFrame}); or a single
#'     list/List.
#'   }
#' }
#'
#' @section Coercion:
#' \describe{
#'   \item{\code{as(x, "RasterLayerList")}, \code{as(x, "PointsLayerList")},
#'     \code{as(x, "ShapesLayerList")}:}{
#'     Coerce \code{list} or \linkS4class{List} to the target class.
#'     For \code{ShapesLayerList}, \code{data.frame} elements are converted to
#'     \code{DataFrame}.
#'   }
#' }
#'
#' @section Displaying:
#' The \code{show()} method prints the class name, length, and for each layer:
#' index, name, element class, and dimensions (raster) or row x col (points/shapes).
#'
#' @note
#' \strong{Terminology alignment with spatialdata:}
#'
#' MASE uses the "Layer" suffix (\code{PointsLayerList}, \code{ShapesLayerList},
#' \code{RasterLayerList}) rather than "Element" to avoid naming collisions with
#' R base classes and Bioconductor conventions. The Python spatialdata package
#' uses "element" as an internal discriminator (e.g., \code{element_type} column,
#' \code{sdata.points}, \code{sdata.shapes}) but does not define List container
#' classes.
#'
#' MASE's \code{element_type} column in \code{spatialMap} corresponds to
#' spatialdata's element type discriminator, routing to the appropriate spatial
#' slot: \code{"points"} → \code{PointsLayerList}, \code{"shapes"} →
#' \code{ShapesLayerList}. See \code{?MultiAssaySpatialExperiment} for details
#' on the \code{spatialMap} linkage model.
#'
#' Other intentional differences: MASE uses \code{spatialPoints()},
#' \code{spatialShapes()}, etc. (with "spatial" prefix) as S4 generic accessors,
#' while spatialdata uses bare attribute access (\code{sdata.points},
#' \code{sdata.shapes}). This follows Bioconductor S4 conventions and avoids
#' conflicts with base R functions.
#'
#' @author Patrick Aboyoun
#'
#' @examples
#' RasterLayerList()
#' RasterLayerList(img = matrix(1:12, 3, 4))
#'
#' pts <- S4Vectors::DataFrame(x = c(1, 2), y = c(3, 4))
#' PointsLayerList(coords = pts)
#'
#' ShapesLayerList()
#' df <- S4Vectors::DataFrame(id = 1:2, value = c("a", "b"))
#' ShapesLayerList(regions = df)
#'
#' @return
#' The constructors \code{RasterLayerList()}, \code{PointsLayerList()}, and
#' \code{ShapesLayerList()} return an object of the corresponding class.
#' Coercion via \code{as()} returns the requested \code{*LayerList}, and the
#' \code{show()} method returns \code{NULL} invisibly.
#'
#' @aliases SpatialLayerList-class
#' @aliases RasterLayerList-class
#' @aliases PointsLayerList-class
#' @aliases ShapesLayerList-class
#'
#' @aliases RasterLayerList
#' @aliases PointsLayerList
#' @aliases ShapesLayerList
#'
#' @aliases show,RasterLayerList-method
#' @aliases show,PointsLayerList-method
#' @aliases show,ShapesLayerList-method
#'
#' @seealso
#' \itemize{
#'   \item \linkS4class{MultiAssaySpatialExperiment} for usage in multi-assay context
#'   \item \linkS4class{SimpleList} for the base class
#' }
#'
#' @keywords classes methods
#'
#' @name SpatialLayerList-class
NULL

#' @importClassesFrom S4Vectors SimpleList
#' @keywords internal
setClass("SpatialLayerList",
    contains = c("SimpleList", "VIRTUAL")
)

#' @exportClass RasterLayerList
setClass("RasterLayerList", contains = "SpatialLayerList")

#' @exportClass PointsLayerList
setClass("PointsLayerList", contains = "SpatialLayerList")

#' @exportClass ShapesLayerList
setClass("ShapesLayerList", contains = "SpatialLayerList")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.checkSpatialLayerListNames <- function(object) {
    if (length(object) == 0L) return(NULL)
    msg <- NULL
    nms <- names(object)
    if (is.null(nms) || any(nms == ""))
        msg <- c(msg, "all elements must be named")
    if (anyDuplicated(nms))
        msg <- c(msg, "element names must be unique")
    msg
}

.validRasterLayerList <- function(object) {
    msg <- .checkSpatialLayerListNames(object)
    if (length(object) == 0L) return(msg)
    for (i in seq_along(object)) {
        el <- object[[i]]
        d <- try(dim(el), silent = TRUE)
        if (inherits(d, "try-error") || is.null(d) || length(d) < 2L)
            msg <- c(msg, sprintf("element [[%s]] must have dim() with length >= 2", i))
    }
    msg
}

#' @importFrom S4Vectors setValidity2
setValidity2("RasterLayerList", .validRasterLayerList)

#' @importClassesFrom S4Vectors DataFrame
.validPointsLayerList <- function(object) {
    msg <- .checkSpatialLayerListNames(object)
    if (length(object) == 0L) return(msg)
    for (i in seq_along(object)) {
        el <- object[[i]]
        if (!is(el, "DataFrame"))
            msg <- c(msg, sprintf("element [[%s]] must be DataFrame", i))
        else if (ncol(el) < 2L)
            msg <- c(msg, sprintf("element [[%s]] must have at least 2 columns", i))
    }
    msg
}

#' @importFrom S4Vectors setValidity2
setValidity2("PointsLayerList", .validPointsLayerList)

#' @importClassesFrom S4Vectors DataFrame
.validShapesLayerList <- function(object) {
    msg <- .checkSpatialLayerListNames(object)
    if (length(object) == 0L) return(msg)
    for (i in seq_along(object)) {
        el <- object[[i]]
        if (!is(el, "DataFrame"))
            msg <- c(msg, sprintf("element [[%s]] must be DataFrame", i))
    }
    msg
}

#' @importFrom S4Vectors setValidity2
setValidity2("ShapesLayerList", .validShapesLayerList)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors
###

#' @importClassesFrom S4Vectors DataFrame List
.normarg_simple_list <- function(...) {
    listData <- list(...)
    if (length(listData) == 1L) {
        if (is.list(listData[[1L]]) && !is.data.frame(listData[[1L]]))
            listData <- listData[[1L]]
        else if (is(listData[[1L]], "List") && !is(listData[[1L]], "DataFrame"))
            listData <- as.list(listData[[1L]])
    }
    if (length(listData) == 0L)
        return(structure(list(), .Names = character()))
    listData
}

#' @importFrom S4Vectors SimpleList
#' @export
RasterLayerList <- function(...) {
    listData <- list(...)
    if (length(listData) == 1L) {
        x <- listData[[1L]]
        if (is(x, "RasterLayerList"))
            return(x)
        if (is(x, "MultiAssaySpatialExperiment") || is(x, "MultiAssayExperiment"))
            stop("MultiAssaySpatialExperiment/MultiAssayExperiment input detected. ",
                 "Did you mean 'spatialImages()' or similar?", call. = FALSE)
    }
    listData <- .normarg_simple_list(...)
    new("RasterLayerList", do.call(SimpleList, listData))
}

#' @importFrom S4Vectors SimpleList
#' @export
PointsLayerList <- function(...) {
    listData <- list(...)
    if (length(listData) == 1L) {
        x <- listData[[1L]]
        if (is(x, "PointsLayerList"))
            return(x)
        if (is(x, "MultiAssaySpatialExperiment") || is(x, "MultiAssayExperiment"))
            stop("MultiAssaySpatialExperiment/MultiAssayExperiment input detected. ",
                 "Did you mean 'spatialPoints()' or similar?", call. = FALSE)
    }
    listData <- .normarg_simple_list(...)
    new("PointsLayerList", do.call(SimpleList, listData))
}

#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom S4Vectors DataFrame SimpleList
#' @export
ShapesLayerList <- function(...) {
    listData <- list(...)
    if (length(listData) == 1L) {
        x <- listData[[1L]]
        if (is(x, "ShapesLayerList"))
            return(x)
        if (is(x, "MultiAssaySpatialExperiment") || is(x, "MultiAssayExperiment"))
            stop("MultiAssaySpatialExperiment/MultiAssayExperiment input detected. ",
                 "Did you mean 'spatialShapes()' or similar?", call. = FALSE)
    }
    listData <- .normarg_simple_list(...)
    for (i in seq_along(listData)) {
        if (is.data.frame(listData[[i]]) && !is(listData[[i]], "DataFrame"))
            listData[[i]] <- DataFrame(listData[[i]])
    }
    new("ShapesLayerList", do.call(SimpleList, listData))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

#' @importFrom S4Vectors SimpleList
setAs("list", "RasterLayerList", function(from) {
    new("RasterLayerList", do.call(SimpleList, from))
})

#' @importFrom S4Vectors SimpleList
setAs("List", "RasterLayerList", function(from) {
    new("RasterLayerList", do.call(SimpleList, as.list(from)))
})

#' @importFrom S4Vectors SimpleList
setAs("list", "PointsLayerList", function(from) {
    new("PointsLayerList", do.call(SimpleList, from))
})

#' @importFrom S4Vectors SimpleList
setAs("List", "PointsLayerList", function(from) {
    new("PointsLayerList", do.call(SimpleList, as.list(from)))
})

#' @importClassesFrom S4Vectors DataFrame
#' @importFrom S4Vectors DataFrame SimpleList
setAs("list", "ShapesLayerList", function(from) {
    from <- lapply(from, function(el) {
        if (is.data.frame(el) && !is(el, "DataFrame"))
            DataFrame(el)
        else el
    })
    new("ShapesLayerList", do.call(SimpleList, from))
})

#' @importClassesFrom S4Vectors DataFrame
#' @importFrom S4Vectors DataFrame SimpleList
setAs("List", "ShapesLayerList", function(from) {
    from <- as.list(from)
    from <- lapply(from, function(el) {
        if (is.data.frame(el) && !is(el, "DataFrame"))
            DataFrame(el)
        else el
    })
    new("ShapesLayerList", do.call(SimpleList, from))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Display
###

.getRasterElDims <- function(el) {
    d <- try(dim(el), silent = TRUE)
    if (inherits(d, "try-error") || is.null(d)) return("?")
    paste(d[seq_len(min(2L, length(d)))], collapse = " x ")
}

#' @export
setMethod("show", "RasterLayerList", function(object) {
    o_class <- class(object)
    o_len <- length(object)
    o_names <- names(object)
    cat(o_class, "of length", o_len, "\n")
    if (o_len > 0L) {
        elem_cl <- vapply(object, function(o) class(o)[[1L]], character(1L))
        elem_dims <- vapply(object, .getRasterElDims, character(1L))
        cat(sprintf("[%i] %s: %s (%s)\n",
            seq_len(o_len), o_names, elem_cl, elem_dims))
    }
})

#' @export
setMethod("show", "PointsLayerList", function(object) {
    o_class <- class(object)
    o_len <- length(object)
    o_names <- names(object)
    cat(o_class, "of length", o_len, "\n")
    if (o_len > 0L) {
        elem_cl <- vapply(object, function(o) class(o)[[1L]], character(1L))
        elem_n <- vapply(object, function(o) paste0(nrow(o), " x ", ncol(o)), character(1L))
        cat(sprintf("[%i] %s: %s (%s)\n",
            seq_len(o_len), o_names, elem_cl, elem_n))
    }
})

#' @export
setMethod("show", "ShapesLayerList", function(object) {
    o_class <- class(object)
    o_len <- length(object)
    o_names <- names(object)
    cat(o_class, "of length", o_len, "\n")
    if (o_len > 0L) {
        elem_cl <- vapply(object, function(o) class(o)[[1L]], character(1L))
        elem_n <- vapply(object, function(o) paste0(nrow(o), " x ", ncol(o)), character(1L))
        cat(sprintf("[%i] %s: %s (%s)\n",
            seq_len(o_len), o_names, elem_cl, elem_n))
    }
})
