### =========================================================================
### spatial-raster-convert.R — label raster ↔ shape conversions
### -------------------------------------------------------------------------

#' @importFrom S4Vectors wmsg
#' @importFrom sf st_as_sf
NULL

.as_label_raster <- function(x) {
    if (!requireNamespace("terra", quietly = TRUE))
        stop(wmsg("Package 'terra' required. ",
                  "Install with: BiocManager::install('terra')"))
    if (is(x, "SpatRaster"))
        return(x)
    if (is.matrix(x) || is.array(x))
        return(terra::rast(x))
    stop(wmsg("'x' must be a SpatRaster, matrix, or array"))
}

.as_shape_table <- function(x) {
    if (is(x, "sf"))
        return(x)
    if (is(x, "data.frame"))
        return(sf::st_as_sf(x))
    stop(wmsg("'x' must be an sf object or data.frame with geometry"))
}

#' Convert label rasters to polygon shapes
#'
#' @description
#' Dissolve contiguous cells with the same label value into polygons. Thin
#' wrapper around \code{\link[terra:as.polygons]{terra::as.polygons}} for
#' segmentation masks stored in \code{spatialLabels()} or vendor TIFFs.
#'
#' @param x A \code{SpatRaster}, 2-D matrix/array, or a single element from a
#'   \linkS4class{RasterLayerList}.
#' @param dissolve Logical; dissolve adjacent cells with the same value
#'   (default \code{TRUE}).
#' @param ... Additional arguments passed to \code{terra::as.polygons}.
#'
#' @return An \pkg{sf} object with dissolved label polygons.
#'
#' @seealso \code{\link{shapesToLabels}}, \code{\link{spatialLabels}}
#'
#' @examples
#' if (requireNamespace("terra", quietly = TRUE)) {
#'     mat <- matrix(c(1L, 1L, 2L, 2L), 2, 2)
#'     labelsToShapes(mat)
#' }
#'
#' @export
labelsToShapes <- function(x, dissolve = TRUE, ...) {
    r <- .as_label_raster(x)
    polys <- terra::as.polygons(r, values = TRUE, dissolve = dissolve, ...)
    sf::st_as_sf(polys)
}

#' Rasterize polygon shapes to a label image
#'
#' @description
#' Burn shape geometries into a label raster. When \code{template} is
#' \code{NULL}, an empty raster is created from the bounding box of \code{x}.
#'
#' @param x An \pkg{sf} object or \code{data.frame} with polygon geometries.
#' @param template Optional \code{SpatRaster} defining extent and resolution.
#' @param field Character; attribute column used as label values (default
#'   \code{"instance_id"}).
#' @param background Numeric background value (default \code{NA}).
#' @param ... Additional arguments passed to \code{terra::rasterize}.
#'
#' @return A \code{SpatRaster} label image.
#'
#' @seealso \code{\link{labelsToShapes}}, \code{\link{spatialShapes}}
#'
#' @examples
#' if (requireNamespace("terra", quietly = TRUE)) {
#'     mat <- matrix(c(1L, 1L, 2L, 2L), 2, 2)
#'     polys <- labelsToShapes(mat)
#'     polys$instance_id <- polys[["lyr.1"]]
#'     shapesToLabels(polys, field = "instance_id", background = 0L)
#' }
#'
#' @export
shapesToLabels <- function(x, template = NULL, field = "instance_id",
                           background = NA, ...) {
    if (!requireNamespace("terra", quietly = TRUE))
        stop(wmsg("Package 'terra' required. ",
                  "Install with: BiocManager::install('terra')"))
    shp <- .as_shape_table(x)
    if (!"geometry" %in% colnames(shp))
        stop(wmsg("'x' must contain a geometry column"))
    if (!field %in% colnames(shp))
        stop(wmsg("Field '", field, "' not found in shape data"))
    v <- terra::vect(shp)
    if (is.null(template)) {
        ext <- terra::ext(v)
        template <- terra::rast(ext, resolution = 1)
    }
    terra::rasterize(v, template, field = field, background = background, ...)
}
