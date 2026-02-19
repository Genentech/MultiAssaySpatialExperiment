### =========================================================================
### annotateWithRegions — point-in-polygon annotation
### -------------------------------------------------------------------------
###
### Add shape instance IDs to spatialMap for points that fall within shapes.
### Enables aggregateByRegion.
###
### -------------------------------------------------------------------------

#' @include MultiAssaySpatialExperiment-class.R
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### annotateWithRegions
###

#' Annotate points with shape regions
#'
#' @description
#' Perform a spatial join between a points layer and a shapes layer. For each
#' point, the matched shape's \code{instance_id} is recorded. This annotation
#' is added to \code{spatialMap} and enables \code{\link{aggregateByRegion}}.
#'
#' @param x A \linkS4class{MultiAssaySpatialExperiment}.
#' @param points Character. Name of the points layer (in \code{spatialPoints(x)}).
#' @param shapes Character. Name of the shapes layer (in \code{spatialShapes(x)}).
#' @param spatialCoordsNames Character vector. Names of the x and y coordinate
#'     columns in the points layer (default \code{c("x", "y")}).
#' @param regionCol Character or \code{NULL}. Name of the column to add to
#'     \code{spatialMap} with the shape instance IDs. Default is the \code{shapes}
#'     name.
#' @param join Function. Spatial join to use, matching \code{\link[sf:st_join]{st_join}}.
#'     Default \code{\link[sf:st_intersects]{st_intersects}} (point-in-polygon).
#'     Use \code{\link[sf:st_nearest_feature]{st_nearest_feature}} to assign each
#'     point to the nearest shape when no shape contains it.
#'
#' @details
#' For each row in the points layer, \code{join} determines which shape (if any)
#' matches. The shape's \code{instance_id} is written to a new column in
#' \code{spatialMap}. Only \code{spatialMap} rows that reference the given
#' \code{points} layer are updated; other rows get \code{NA} for the new column.
#'
#' With \code{join = st_intersects} (default), points not falling within any
#' shape receive \code{NA}. With \code{join = st_nearest_feature}, every point
#' is assigned to its nearest shape. The user is responsible for ensuring
#' \code{instance_id} types match between points and shapes (no coercion).
#'
#' @section Polymorphism:
#' Uses S3/S4 generics only. Table-like elements (\code{spatialMap}, points, shapes)
#' use \code{nrow}, \code{colnames}, \code{[[}, and \code{[}. Spatial operations use
#' the supplied \code{join} function, \code{\link[sf:st_as_sfc]{st_as_sfc}}, and
#' \code{\link[sf:st_crs]{st_crs}}. Layer lists use \code{names}. Consistent with
#' \code{\link{subsetByPolygon}}; see the Subset vignette.
#'
#' @return
#' \code{x} with \code{spatialMap} updated to include the new annotation column.
#' If \code{spatialMap} is \code{NULL} or the required layers are missing,
#' \code{x} is returned unchanged.
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
#'         dimnames = list(NULL, c("A", "B", "C")))),
#'     colData = S4Vectors::DataFrame(row.names = "s1"),
#'     sampleMap = S4Vectors::DataFrame(
#'         assay = "assay1", primary = "s1", colname = c("A", "B", "C")),
#'     points = PointsLayerList(centroids = pts),
#'     shapes = ShapesLayerList(cells = shp_df),
#'     spatialMap = S4Vectors::DataFrame(
#'         assay = "assay1", colname = c("A", "B", "C"),
#'         region = "centroids", instance_id = c("A", "B", "C")))
#' mase <- annotateWithRegions(mase, points = "centroids", shapes = "cells")
#' spatialMap(mase)
#'
#' @aliases
#' annotateWithRegions
#' annotateWithRegions,MultiAssaySpatialExperiment-method
#'
#' @export
#' @importFrom sf st_intersects st_nearest_feature
setGeneric("annotateWithRegions",
    function(x, points, shapes, spatialCoordsNames = c("x", "y"), regionCol = NULL,
        join = st_intersects)
        standardGeneric("annotateWithRegions"))

#' @importFrom sf st_as_sfc st_crs st_intersects
.pointToShapeMapping <- function(pt, shp, x_col, y_col, geom_col = "geometry",
                                 id_col = "instance_id", join = st_intersects) {
    if (is.null(pt) || nrow(pt) == 0L || is.null(shp) || nrow(shp) == 0L)
        return(structure(logical(0), names = character(0)))
    if (!x_col %in% colnames(pt) || !y_col %in% colnames(pt))
        return(structure(logical(0), names = character(0)))
    if (!geom_col %in% colnames(shp) || !id_col %in% colnames(shp))
        return(structure(logical(0), names = character(0)))
    pts_sfc <- st_as_sfc(
        paste0("POINT(", pt[[x_col]], " ", pt[[y_col]], ")"),
        crs = st_crs(shp[[geom_col]]))
    shp_geom <- shp[[geom_col]]
    res <- join(pts_sfc, shp_geom)
    n_pts <- nrow(pt)
    shp_ids <- shp[[id_col]]
    pt_to_shp <- vector("list", n_pts)
    if (inherits(res, "sgbp") || (is.list(res) && length(res) == n_pts)) {
        for (i in seq_len(n_pts)) {
            idx <- res[[i]]
            pt_to_shp[[i]] <- if (length(idx) > 0L) shp_ids[idx[1L]] else NA
        }
    } else if (is.integer(res) || is.numeric(res)) {
        for (i in seq_len(n_pts)) {
            idx <- res[i]
            if (!is.na(idx) && idx >= 1L && idx <= length(shp_ids))
                pt_to_shp[[i]] <- shp_ids[idx]
            else
                pt_to_shp[[i]] <- NA
        }
    } else {
        stop("'join' must return an sgbp list or an integer vector")
    }
    pt_to_shp <- unlist(pt_to_shp, use.names = FALSE)
    pt_id <- if (id_col %in% colnames(pt)) pt[[id_col]] else seq_len(n_pts)
    structure(pt_to_shp, names = as.character(pt_id))
}

#' @importFrom MultiAssayExperiment experiments
#' @exportMethod annotateWithRegions
setMethod("annotateWithRegions", "MultiAssaySpatialExperiment",
function(x, points, shapes, spatialCoordsNames = c("x", "y"), regionCol = NULL,
    join = st_intersects) {
    pts_el <- spatialPoints(x)[[points]]
    shps_el <- spatialShapes(x)[[shapes]]
    spmap <- spatialMap(x)

    if (is.null(pts_el) || is.null(shps_el) || is.null(spmap) || nrow(spmap) == 0L)
        return(x)
    if (!points %in% names(spatialPoints(x)) || !shapes %in% names(spatialShapes(x)))
        return(x)
    if (length(spatialCoordsNames) < 2L)
        return(x)
    x_col <- spatialCoordsNames[1L]
    y_col <- spatialCoordsNames[2L]

    if (is.null(regionCol))
        regionCol <- shapes

    pt_to_shp <- .pointToShapeMapping(pts_el, shps_el, x_col, y_col, join = join)
    if (length(pt_to_shp) == 0L)
        return(x)

    keep <- spmap[["region"]] == points
    if (!any(keep))
        return(x)

    annot <- pt_to_shp[as.character(spmap[["instance_id"]])]
    new_col <- spmap[[regionCol]]
    if (is.null(new_col))
        new_col <- rep(NA, nrow(spmap))
    new_col[keep] <- annot[keep]
    spmap[[regionCol]] <- new_col

    replaceSlots(x, spatialMap = spmap, check = FALSE)
})
