### =========================================================================
### spatial-operations.R — S4 generics for spatial operations on DataFrame
### -------------------------------------------------------------------------
###
### Exported S4 generics following the S4Vectors pattern.  MASE defines
### these generics and provides DataFrame methods (sf path).  BiocDuckDB
### imports the generics and provides DuckDBDataFrame methods (SQL path).
###
### spatialOverlaps  — predicate (like is.na, duplicated); returns logical
### spatialMatch     — relational (like match); returns integer positions
###
### =========================================================================

#' @include SpatialGenerics.R
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### spatialOverlaps
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Spatial overlap predicate for DataFrame
#'
#' @description
#' Test which rows of a \code{\linkS4class{DataFrame}} spatially intersect a
#' region.  Analogous to \code{is.na} or \code{duplicated} from
#' \pkg{S4Vectors}: returns a logical vector that can be used for subsetting
#' via \code{x[spatialOverlaps(x, y, ...), ]}.
#'
#' @param x A \code{\linkS4class{DataFrame}} with spatial data (coordinate
#'   columns or a geometry column).
#' @param y A spatial region to test against.  Typically an \code{sfc},
#'   \code{sfg}, or similar geometry object from \pkg{sf}.
#' @param ... Additional arguments passed to methods.  Common arguments:
#'   \describe{
#'     \item{\code{coords}}{Character vector of length 2 naming the
#'       coordinate columns in \code{x} (e.g., \code{c("x", "y")}).
#'       When provided, point geometries are constructed from these columns
#'       and tested for intersection with \code{y}.}
#'     \item{\code{geom}}{Character; name of the geometry column in
#'       \code{x} (default \code{"geometry"}).  Used when \code{coords}
#'       is not provided.}
#'   }
#'
#' @return A logical vector of length \code{nrow(x)}.  \code{TRUE} for rows
#'   whose spatial data intersects \code{y}.
#'
#' @details
#' The \code{DataFrame} method uses \pkg{sf} functions (\code{st_as_sfc},
#' \code{st_intersects}).  Downstream packages (e.g., BiocDuckDB) can define
#' methods for their own \code{DataFrame} subclasses to use alternative
#' backends (e.g., DuckDB SQL).
#'
#' Two modes are supported based on the arguments:
#' \describe{
#'   \item{Coordinate mode}{When \code{coords} is provided, point geometries
#'     are constructed from those columns and tested against \code{y}.}
#'   \item{Geometry mode}{When \code{coords} is \code{NULL} (the default),
#'     the column named by \code{sf_column_name} is used directly.}
#' }
#'
#' @seealso \code{\link{spatialMatch}} for the relational analogue (like
#'   \code{match}).
#'
#' @author Patrick Aboyoun
#'
#' @examples
#' library(sf)
#' pts <- S4Vectors::DataFrame(x = c(1, 5, 10), y = c(1, 5, 10))
#' polygon <- st_as_sfc("POLYGON((0 0, 6 0, 6 6, 0 6, 0 0))")
#' spatialOverlaps(pts, polygon, coords = c("x", "y"))
#' ## TRUE  TRUE FALSE
#'
#' @aliases
#' spatialOverlaps
#' spatialOverlaps,DataFrame-method
#'
#' @name spatialOverlaps
NULL

#' @importClassesFrom S4Vectors DataFrame
#' @export
setGeneric("spatialOverlaps", signature = "x",
    function(x, y, ...) standardGeneric("spatialOverlaps"))

#' @importFrom sf st_as_sfc st_intersects
#' @exportMethod spatialOverlaps
setMethod("spatialOverlaps", "DataFrame",
    function(x, y, coords = NULL, geom = "geometry", ...) {
        if (!is.null(coords)) {
            pts_sfc <- st_as_sfc(
                paste0("POINT(", x[[coords[1L]]], " ", x[[coords[2L]]], ")"))
            lengths(st_intersects(pts_sfc, y)) > 0L
        } else {
            geom_col <- x[[geom]]
            lengths(st_intersects(geom_col, y)) > 0L
        }
    })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### spatialMatch
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Spatial match for DataFrame
#'
#' @description
#' For each row of \code{x}, find the first matching row in \code{table}
#' based on a spatial relationship.  Analogous to \code{match(x, table)}
#' from \pkg{S4Vectors}: returns an integer vector of positions in
#' \code{table}.
#'
#' @param x A \code{\linkS4class{DataFrame}} with coordinate columns
#'   representing point locations.
#' @param table A \code{\linkS4class{DataFrame}} with a geometry column
#'   representing spatial regions.
#' @param ... Additional arguments passed to methods.  Common arguments:
#'   \describe{
#'     \item{\code{coords}}{Character vector of length 2 naming the
#'       coordinate columns in \code{x}.}
#'     \item{\code{geom}}{Character; name of the geometry column in
#'       \code{table} (default \code{"geometry"}).}
#'     \item{\code{join}}{A spatial join function (default
#'       \code{sf::st_intersects}).  Must accept two geometry arguments and
#'       return either an sgbp list or an integer vector.}
#'   }
#'
#' @return An integer vector of length \code{nrow(x)}.  Each element is the
#'   row position in \code{table} of the first spatial match, or \code{NA}
#'   if no match is found.
#'
#' @details
#' The \code{DataFrame} method builds point geometries from the coordinate
#' columns of \code{x}, applies the \code{join} function against the
#' geometry column of \code{table}, and returns the first match per row.
#'
#' Downstream packages (e.g., BiocDuckDB) can define methods for their own
#' \code{DataFrame} subclasses to push the spatial join to an alternative
#' backend (e.g., a DuckDB SQL spatial join).
#'
#' @seealso \code{\link{spatialOverlaps}} for the predicate analogue (like
#'   \code{is.na}).
#'
#' @author Patrick Aboyoun
#'
#' @examples
#' library(sf)
#' pts <- S4Vectors::DataFrame(x = c(1, 5, 10), y = c(1, 5, 10))
#' shapes <- S4Vectors::DataFrame(
#'     geometry = st_as_sfc(c(
#'         "POLYGON((0 0, 6 0, 6 6, 0 6, 0 0))",
#'         "POLYGON((7 7, 12 7, 12 12, 7 12, 7 7))")))
#' spatialMatch(pts, shapes, coords = c("x", "y"))
#' ## 1  1  2
#'
#' @aliases
#' spatialMatch
#' spatialMatch,DataFrame,DataFrame-method
#'
#' @name spatialMatch
NULL

#' @importClassesFrom S4Vectors DataFrame
#' @export
setGeneric("spatialMatch", signature = c("x", "table"),
    function(x, table, ...) standardGeneric("spatialMatch"))

#' @importFrom sf st_as_sfc st_intersects
#' @exportMethod spatialMatch
setMethod("spatialMatch", c("DataFrame", "DataFrame"),
    function(x, table, coords, geom = "geometry",
             join = NULL, ...) {
        if (is.null(join))
            join <- st_intersects
        pts_sfc <- st_as_sfc(
            paste0("POINT(", x[[coords[1L]]], " ", x[[coords[2L]]], ")"))
        tbl_geom <- table[[geom]]
        res <- join(pts_sfc, tbl_geom)
        n <- nrow(x)
        out <- integer(n)
        if (inherits(res, "sgbp") || (is.list(res) && length(res) == n)) {
            for (i in seq_len(n)) {
                idx <- res[[i]]
                out[i] <- if (length(idx) > 0L) idx[1L] else NA_integer_
            }
        } else if (is.integer(res) || is.numeric(res)) {
            for (i in seq_len(n)) {
                idx <- res[i]
                out[i] <- if (!is.na(idx) && idx >= 1L && idx <= nrow(table))
                    as.integer(idx) else NA_integer_
            }
        } else {
            stop("'join' must return an sgbp list or an integer vector")
        }
        out
    })
