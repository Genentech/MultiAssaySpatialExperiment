### =========================================================================
### MultiAssaySpatialExperiment-helpers.R — prepMultiAssay wrapper for MASE
### -------------------------------------------------------------------------

#' @include MultiAssaySpatialExperiment-spatialMap.R
#' @include MultiAssaySpatialExperiment-class.R
NULL

#' Prepare components for a MultiAssaySpatialExperiment
#'
#' @description
#' Wraps \code{\link[MultiAssayExperiment]{prepMultiAssay}} and harmonizes spatial
#' slots. Returns a \code{list} suitable for
#' \code{\link{MultiAssaySpatialExperiment}} or \code{do.call}.
#'
#' @param ExperimentList A \link[MultiAssayExperiment]{ExperimentList} or named
#'   list of assays passed to \code{\link[MultiAssayExperiment]{prepMultiAssay}}.
#' @param colData A \linkS4class{DataFrame} of specimen-level metadata (one row
#'   per \code{sampleMap$primary}), passed to \code{prepMultiAssay}.
#' @param sampleMap A \linkS4class{DataFrame} with MAE columns \code{assay},
#'   \code{primary}, and \code{colname}, passed to \code{prepMultiAssay}.
#' @param points A \linkS4class{PointsLayerList} or named list of point layers.
#' @param shapes A \linkS4class{ShapesLayerList} or named list of shape layers.
#' @param images A \linkS4class{RasterLayerList} of image layers (default empty).
#' @param labels A \linkS4class{RasterLayerList} of label layers (default empty).
#' @param imgData Optional specimen-level image metadata \linkS4class{DataFrame}.
#' @param spatialMap Optional instance-level mapping \linkS4class{DataFrame}.
#' @param ... Additional arguments passed to \code{prepMultiAssay} (\code{metadata},
#'   \code{drops}).
#'
#' @return A \code{list} with elements \code{experiments}, \code{colData},
#'   \code{sampleMap}, \code{metadata}, \code{drops}, and spatial slot components.
#'
#' @seealso \code{\link{buildSpatialMap}},
#'   \code{\link[MultiAssayExperiment]{prepMultiAssay}}
#'
#' @export
#' @importFrom MultiAssayExperiment prepMultiAssay
prepMASE <-
function(ExperimentList, colData, sampleMap, points = PointsLayerList(),
         shapes = ShapesLayerList(), images = RasterLayerList(),
         labels = RasterLayerList(), imgData = NULL, spatialMap = NULL, ...)
{
    prepared <- prepMultiAssay(ExperimentList, colData, sampleMap, ...)
    meta <- prepared$metadata
    drops <- if (is.list(meta) && !is.null(meta[["drops"]])) meta[["drops"]] else list()
    if (!is.list(drops)) drops <- list()
    pts <- if (is(points, "PointsLayerList")) points else PointsLayerList(points)
    shps <- if (is(shapes, "ShapesLayerList")) shapes else ShapesLayerList(shapes)
    if (!is.null(spatialMap) && nrow(spatialMap) > 0L) {
        harm <- .harmonizeSpatialMap(spatialMap, prepared$sampleMap, pts, shps)
        spatialMap <- harm$spatialMap
        if (length(harm$drops))
            drops <- c(drops, harm$drops)
    }
    prepared$points <- pts
    prepared$shapes <- shps
    prepared$images <- if (is(images, "RasterLayerList")) images else
        RasterLayerList(images)
    prepared$labels <- if (is(labels, "RasterLayerList")) labels else
        RasterLayerList(labels)
    prepared$imgData <- imgData
    prepared$spatialMap <- spatialMap
    if (is.list(prepared$metadata))
        prepared$metadata[["drops"]] <- drops
    prepared
}
