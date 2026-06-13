### =========================================================================
### MultiAssaySpatialExperiment - SpatialExperiment integration
### -------------------------------------------------------------------------
###
### Coercion to/from SpatialExperiment; methods for spatialCoords, getImg,
### imgRaster, scaleFactors, addImg, rmvImg, etc.
###
### -------------------------------------------------------------------------

#' @include MultiAssaySpatialExperiment-class.R
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion and SpatialExperiment methods
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' MultiAssaySpatialExperiment - SpatialExperiment integration
#'
#' @description
#' Coercion between \linkS4class{SpatialExperiment} and
#' \linkS4class{MultiAssaySpatialExperiment}, plus methods for
#' \code{spatialCoords}, \code{imgData}, \code{scaleFactors}, \code{getImg},
#' \code{addImg}, \code{rmvImg}, and related SpatialExperiment-derived
#' functionality.
#'
#' @section Coercion:
#' \describe{
#'   \item{\code{as(spe, "MultiAssaySpatialExperiment")}:}{
#'     Wraps a \linkS4class{SpatialExperiment} in MASE with a single assay.
#'     \code{spatialCoords} are copied to points; \code{imgData} is passed
#'     through.
#'   }
#'   \item{\code{as(mase, "SpatialExperiment")}:}{
#'     Requires exactly one assay that is or can be coerced to
#'     \linkS4class{SpatialExperiment}.
#'   }
#' }
#'
#' @section spatialCoords and scaleFactors:
#' \describe{
#'   \item{\code{spatialCoords(x, assay = 1L)}:}{
#'     Returns spatial coordinates from the assay (if a SpatialExperiment) or
#'     from \code{spatialPoints(x)}.
#'   }
#'   \item{\code{spatialCoordsNames(x)}:}{
#'     Returns coordinate column names from the first assay or points layer.
#'   }
#'   \item{\code{scaleFactors(x, sample_id = TRUE, image_id = TRUE)}:}{
#'     Returns scale factors from \code{imgData}.
#'   }
#' }
#'
#' @section imgData methods:
#' \describe{
#'   \item{\code{getImg(x, sample_id = NULL, image_id = NULL)}:}{
#'     Retrieves image(s) from \code{imgData}.
#'   }
#'   \item{\code{addImg(x, imageSource, scaleFactor, sample_id, image_id,
#'     load = TRUE)}:}{
#'     Adds an image to \code{imgData}.
#'   }
#'   \item{\code{rmvImg(x, sample_id = NULL, image_id = NULL)}:}{
#'     Removes image(s) from \code{imgData}.
#'   }
#'   \item{\code{imgSource(x, sample_id = NULL, image_id = NULL, path = FALSE)}:}{
#'     Returns image source path(s).
#'   }
#'   \item{\code{imgRaster(x, sample_id = NULL, image_id = NULL)}:}{
#'     Returns raster image(s).
#'   }
#'   \item{\code{rotateImg(x, sample_id = NULL, image_id = NULL, degrees = 90)}:}{
#'     Rotates image(s) in place.
#'   }
#'   \item{\code{mirrorImg(x, sample_id = NULL, image_id = NULL,
#'     axis = c("h", "v"))}:}{
#'     Mirrors image(s) in place.
#'   }
#' }
#'
#' @section molecules:
#' \describe{
#'   \item{\code{molecules(x, ...)}:}{
#'     Returns molecule data from the first assay if it is a SpatialExperiment
#'     with a \code{molecules} assay; otherwise \code{NULL}.
#'   }
#' }
#'
#' @author Patrick Aboyoun
#'
#' @aliases coerce,SpatialExperiment,MultiAssaySpatialExperiment-method
#' @aliases coerce,MultiAssaySpatialExperiment,SpatialExperiment-method
#'
#' @aliases spatialCoords,MultiAssaySpatialExperiment-method
#' @aliases spatialCoordsNames,MultiAssaySpatialExperiment-method
#' @aliases scaleFactors,MultiAssaySpatialExperiment-method
#'
#' @aliases getImg,MultiAssaySpatialExperiment-method
#' @aliases addImg,MultiAssaySpatialExperiment-method
#' @aliases rmvImg,MultiAssaySpatialExperiment-method
#' @aliases imgSource,MultiAssaySpatialExperiment-method
#' @aliases imgRaster,MultiAssaySpatialExperiment-method
#' @aliases rotateImg,MultiAssaySpatialExperiment-method
#' @aliases mirrorImg,MultiAssaySpatialExperiment-method
#'
#' @aliases molecules,MultiAssaySpatialExperiment-method
#'
#' @seealso
#' \itemize{
#'   \item \linkS4class{MultiAssaySpatialExperiment} for the main class
#'   \item \linkS4class{SpatialExperiment} for single-assay spatial data
#' }
#'
#' @keywords methods
#'
#' @name MultiAssaySpatialExperiment-SpatialExperiment
NULL

#' @importClassesFrom SpatialExperiment SpatialExperiment
#' @importFrom MultiAssayExperiment ExperimentList MultiAssayExperiment
#' @importFrom SpatialExperiment imgData spatialCoords
#' @importFrom S4Vectors DataFrame
#' @importFrom stats setNames
#' @export
setAs("SpatialExperiment", "MultiAssaySpatialExperiment", function(from) {
    assay_name <- "spatial"
    n <- ncol(from)
    primaries <- unique(from$sample_id)
    if (length(primaries) == 0L)
        primaries <- "sample01"
    cd <- DataFrame(row.names = primaries)
    colnms <- colnames(from)
    if (is.null(colnms) || length(colnms) != n)
        colnms <- if (n > 0L) paste0("cell", seq_len(n)) else character(0L)
    if (is.null(colnames(from)) || !identical(colnames(from), colnms)) {
        from <- initialize(from)
        colnames(from) <- colnms
    }
    experiments <- ExperimentList(setNames(list(from), assay_name))
    sids <- from$sample_id
    if (is.null(sids) || length(sids) != n)
        sids <- rep(primaries[1L], n)
    sm <- DataFrame(assay = factor(rep(assay_name, n), assay_name),
                    primary = as.character(sids),
                    colname = colnms)
    mae <- MultiAssayExperiment(experiments = experiments, colData = cd,
                                sampleMap = sm)
    pts <- PointsLayerList()
    if (ncol(from) > 0L) {
        scoords <- spatialCoords(from)
        if (!is.null(scoords) && nrow(scoords) > 0L) {
            pt_df <- DataFrame(as.matrix(scoords))
            rn <- rownames(scoords)
            pt_df[["instance_id"]] <- if (!is.null(rn) && length(rn) > 0L) rn else colnms
            pts <- PointsLayerList(coordinates = pt_df)
        }
    }
    img <- imgData(from)
    new("MultiAssaySpatialExperiment",
        mae,
        images = RasterLayerList(),
        labels = RasterLayerList(),
        points = pts,
        shapes = ShapesLayerList(),
        imgData = img,
        spatialMap = NULL)
})

#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom SpatialExperiment SpatialExperiment
#' @importFrom MultiAssayExperiment experiments
#' @importFrom SpatialExperiment imgData imgData<-
#' @export
setAs("MultiAssaySpatialExperiment", "SpatialExperiment", function(from) {
    exps <- experiments(from)
    if (length(exps) != 1L)
        stop("Coercion to SpatialExperiment requires exactly one assay; ",
             "found ", length(exps))
    exp <- exps[[1L]]
    if (is(exp, "SpatialExperiment")) {
        spe <- exp
    } else if (is(exp, "SingleCellExperiment")) {
        spe <- as(exp, "SpatialExperiment")
    } else {
        stop("The single assay must be a SpatialExperiment or ",
             "SingleCellExperiment; found ", class(exp)[1L])
    }
    mase_img <- imgData(from)
    if (!is.null(mase_img))
        imgData(spe) <- mase_img
    spe
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper: imgData row index
###

#' @importFrom SpatialExperiment imgData
.get_img_idx_mase <- function(x, sample_id = NULL, image_id = NULL) {
    img <- imgData(x)
    if (is.null(img) || nrow(img) == 0L)
        stop("'imgData' is empty")
    for (i in c("sample_id", "image_id")) {
        j <- get(i)
        if (is.factor(j) || is.numeric(j))
            assign(i, as.character(j))
        if (!(is.null(j) || j %in% img[[i]] ||
            (length(j) == 1L && is.logical(j))))
            stop(sprintf(
                "'%s' invalid; should be NULL, TRUE/FALSE, or matching entries in imgData(.)$%s",
                i, i))
    }
    if (is.character(sample_id) && is.character(image_id)) {
        sid <- img$sample_id == sample_id
        iid <- img$image_id == image_id
    } else if (isTRUE(sample_id) && isTRUE(image_id)) {
        sid <- iid <- !logical(nrow(img))
    } else if (is.null(sample_id) && is.null(image_id)) {
        sid <- iid <- diag(nrow(img))[1L, ]
    } else if (is.character(sample_id) && isTRUE(image_id)) {
        sid <- img$sample_id == sample_id
        iid <- !logical(nrow(img))
    } else if (is.character(image_id) && isTRUE(sample_id)) {
        iid <- img$image_id == image_id
        sid <- !logical(nrow(img))
    } else if (is.character(sample_id) && is.null(image_id)) {
        sid <- img$sample_id == sample_id
        iid <- diag(nrow(img))[which(sid)[1L], ]
    } else if (is.character(image_id) && is.null(sample_id)) {
        iid <- img$image_id == image_id
        sid <- diag(nrow(img))[which(iid)[1L], ]
    } else if (isTRUE(sample_id) && is.null(image_id)) {
        iid <- match(unique(img$sample_id), img$sample_id)
        iid <- colSums(diag(nrow(img))[iid, , drop = FALSE])
        sid <- !logical(nrow(img))
    } else if (isTRUE(image_id) && is.null(sample_id)) {
        sid <- match(unique(img$image_id), img$image_id)
        sid <- colSums(diag(nrow(img))[sid, , drop = FALSE])
        iid <- !logical(nrow(img))
    } else {
        stop("invalid 'sample_id' / 'image_id' combination")
    }
    if (!any(idx <- sid & iid))
        stop("No 'imgData' entry matched the specified 'image_id' and 'sample_id'")
    which(idx)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### spatialCoords and scaleFactors
###

#' @importClassesFrom SpatialExperiment SpatialExperiment
#' @importFrom MultiAssayExperiment experiments
#' @importFrom SpatialExperiment spatialCoords
#' @export
setMethod("spatialCoords", "MultiAssaySpatialExperiment",
    function(x, assay = 1L) {
        exps <- experiments(x)
        if (length(exps) == 0L)
            return(matrix(numeric(0L), 0L, 0L))
        exp <- exps[[assay]]
        if (is(exp, "SpatialExperiment"))
            return(spatialCoords(exp))
        pts <- spatialPoints(x)
        if (length(pts) > 0L && "coordinates" %in% spatialPointNames(x)) {
            coords <- pts[["coordinates"]]
            coord_cols <- setdiff(colnames(coords), "instance_id")
            if (length(coord_cols) >= 2L)
                return(as.matrix(coords[, coord_cols]))
        }
        matrix(numeric(0L), ncol(exp), 0L)
    }
)

#' @importClassesFrom SpatialExperiment SpatialExperiment
#' @importFrom MultiAssayExperiment experiments
#' @importFrom SpatialExperiment spatialCoordsNames
#' @export
setMethod("spatialCoordsNames", "MultiAssaySpatialExperiment",
    function(x) {
        exps <- experiments(x)
        if (length(exps) == 0L)
            return(character(0L))
        exp <- exps[[1L]]
        if (is(exp, "SpatialExperiment")) {
            nms <- spatialCoordsNames(exp)
            if (length(nms) > 0L)
                return(nms)
        }
        pts <- spatialPoints(x)
        if (length(pts) > 0L && "coordinates" %in% spatialPointNames(x)) {
            coord_cols <- setdiff(colnames(pts[["coordinates"]]), "instance_id")
            if (length(coord_cols) >= 2L)
                return(coord_cols)
        }
        c("x", "y")
    }
)

#' @importFrom SpatialExperiment imgData scaleFactors
#' @export
setMethod("scaleFactors", "MultiAssaySpatialExperiment",
    function(x, sample_id = TRUE, image_id = TRUE) {
        if (is.null(imgData(x)))
            stop("'imgData' is NULL")
        idx <- .get_img_idx_mase(x, sample_id, image_id)
        imgData(x)$scaleFactor[idx]
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### imgData methods
###

#' @importFrom S4Vectors isEmpty
#' @importFrom SpatialExperiment getImg imgData
#' @export
setMethod("getImg", "MultiAssaySpatialExperiment",
    function(x, sample_id = NULL, image_id = NULL) {
        if (isEmpty(imgData(x)))
            return(NULL)
        spi <- imgData(x)$data
        idx <- .get_img_idx_mase(x, sample_id, image_id)
        if (length(idx) == 1L) spi[[idx]] else spi[idx]
    }
)

#' @importFrom SpatialExperiment addImg imgData imgData<-
#' @export
setMethod("addImg", "MultiAssaySpatialExperiment",
    function(x, imageSource, scaleFactor, sample_id, image_id, load = TRUE) {
        stopifnot(
            is.numeric(scaleFactor), length(scaleFactor) == 1L,
            is.character(sample_id), length(sample_id) == 1L,
            sample_id %in% unique(sampleMap(x)$primary),
            is.character(image_id), length(image_id) == 1L,
            is.logical(load), length(load) == 1L)
        idx <- tryCatch(
            error = function(e) e,
            .get_img_idx_mase(x, sample_id, image_id))
        if (!inherits(idx, "error"))
            stop("'imgData' already contains an entry with ",
                sprintf("'image_id = %s' and 'sample_id = %s'",
                    dQuote(image_id), dQuote(sample_id)))
        df <- SpatialExperiment:::.get_imgData(
            imageSource, scaleFactor, sample_id, image_id, load)
        idf <- imgData(x)
        idf <- if (is.null(idf) || nrow(idf) == 0L) df else rbind(idf, df)
        imgData(x) <- idf
        x
    }
)

#' @importFrom SpatialExperiment imgData imgData<- rmvImg
#' @export
setMethod("rmvImg", "MultiAssaySpatialExperiment",
    function(x, sample_id = NULL, image_id = NULL) {
        idf <- imgData(x)
        if (is.null(idf) || nrow(idf) == 0L)
            return(x)
        idx <- .get_img_idx_mase(x, sample_id, image_id)
        imgData(x) <- idf[-idx, , drop = FALSE]
        x
    }
)

#' @importFrom SpatialExperiment getImg imgSource
#' @export
setMethod("imgSource", "MultiAssaySpatialExperiment",
    function(x, sample_id = NULL, image_id = NULL, path = FALSE) {
        spi <- getImg(x, sample_id, image_id)
        if (is.null(spi))
            return(NULL)
        if (!is.list(spi))
            spi <- list(spi)
        vapply(spi, function(.x) imgSource(.x, path), character(1L))
    }
)

#' @importFrom SpatialExperiment getImg imgRaster
#' @export
setMethod("imgRaster", "MultiAssaySpatialExperiment",
    function(x, sample_id = NULL, image_id = NULL) {
        spi <- getImg(x, sample_id, image_id)
        if (is.null(spi))
            return(NULL)
        if (is.list(spi))
            lapply(spi, imgRaster)
        else
            imgRaster(spi)
    }
)

#' @importFrom SpatialExperiment getImg imgData imgData<- rotateImg
#' @export
setMethod("rotateImg", "MultiAssaySpatialExperiment",
    function(x, sample_id = NULL, image_id = NULL, degrees = 90) {
        old <- getImg(x, sample_id, image_id)
        if (!is.null(old)) {
            if (!is.list(old))
                old <- list(old)
            new <- lapply(old, rotateImg, degrees = degrees)
            idx <- .get_img_idx_mase(x, sample_id, image_id)
            idf <- imgData(x)
            idf$data[idx] <- new
            imgData(x) <- idf
        }
        x
    }
)

#' @importFrom SpatialExperiment getImg imgData imgData<- mirrorImg
#' @export
setMethod("mirrorImg", "MultiAssaySpatialExperiment",
    function(x, sample_id = NULL, image_id = NULL, axis = c("h", "v")) {
        axis <- match.arg(axis)
        old <- getImg(x, sample_id, image_id)
        if (!is.null(old)) {
            if (!is.list(old))
                old <- list(old)
            new <- lapply(old, mirrorImg, axis = axis)
            idx <- .get_img_idx_mase(x, sample_id, image_id)
            idf <- imgData(x)
            idf$data[idx] <- new
            imgData(x) <- idf
        }
        x
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### molecules
###

#' @importClassesFrom SpatialExperiment SpatialExperiment
#' @importFrom MultiAssayExperiment experiments
#' @importFrom SpatialExperiment molecules
#' @importFrom SummarizedExperiment assayNames
#' @export
setMethod("molecules", "MultiAssaySpatialExperiment",
    function(x, ...) {
        exps <- experiments(x)
        if (length(exps) == 0L)
            return(NULL)
        exp <- exps[[1L]]
        if (is(exp, "SpatialExperiment") &&
            "molecules" %in% assayNames(exp))
            molecules(exp, ...)
        else
            NULL
    }
)
