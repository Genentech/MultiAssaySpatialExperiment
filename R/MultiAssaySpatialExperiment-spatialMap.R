### =========================================================================
### MultiAssaySpatialExperiment-spatialMap.R — spatialMap FK utilities
### -------------------------------------------------------------------------
### Shared constants, construction (buildSpatialMap), harmonization, and checks.
### -------------------------------------------------------------------------

#' @include SpatialLayerList-class.R
NULL

SPATIAL_MAP_COLS <- c("assay", "colname", "element_type", "region", "instance_id")
VALID_ELEMENT_TYPES <- c("points", "shapes")

.filterSpatialMapByInstances <- function(spatialMap, instance_ids) {
    if (is.null(spatialMap) || nrow(spatialMap) == 0L)
        return(spatialMap)
    instance_ids <- as.character(instance_ids)
    spatialMap[as.character(spatialMap[["instance_id"]]) %in% instance_ids, ,
               drop = FALSE]
}

.rbindSpatialMaps <- function(...) {
    maps <- list(...)
    maps <- Filter(function(m) !is.null(m) && nrow(m) > 0L, maps)
    if (!length(maps))
        return(NULL)
    do.call(rbind, maps)
}

.transcriptSpatialMapRow <- function(assay, n, region = "transcripts") {
    assay_chr <- as.character(assay)[[1L]]
    DataFrame(assay = factor(rep(assay_chr, n), assay_chr),
              colname = NA_character_,
              element_type = "points",
              region = rep(region, n),
              instance_id = seq_len(n))
}

.harmonizeSpatialMap <- function(spatialMap, sampleMap, points, shapes) {
    drops <- list()
    if (is.null(spatialMap) || nrow(spatialMap) == 0L)
        return(list(spatialMap = spatialMap, drops = drops))
    if (!all(SPATIAL_MAP_COLS %in% colnames(spatialMap)))
        stop(wmsg("'spatialMap' must have columns: ",
                  paste(SPATIAL_MAP_COLS, collapse = ", ")))
    n_before <- nrow(spatialMap)
    if (nrow(sampleMap) > 0L &&
            all(c("assay", "colname") %in% colnames(sampleMap))) {
        sm_key <- paste(sampleMap[["assay"]], sampleMap[["colname"]], sep = "\r")
        linked <- !is.na(spatialMap[["colname"]])
        sp_key <- paste(spatialMap[["assay"]], spatialMap[["colname"]], sep = "\r")
        keep <- !linked | sp_key %in% sm_key
        if (!all(keep)) {
            drops[["spatialMap_sampleMap"]] <- spatialMap[!keep, , drop = FALSE]
            spatialMap <- spatialMap[keep, , drop = FALSE]
        }
    }
    if (nrow(spatialMap) == 0L)
        return(list(spatialMap = spatialMap, drops = drops))
    orphan_inst <- logical(nrow(spatialMap))
    for (i in seq_len(nrow(spatialMap))) {
        et <- spatialMap[["element_type"]][[i]]
        reg <- spatialMap[["region"]][[i]]
        inst <- spatialMap[["instance_id"]][[i]]
        if (!et %in% VALID_ELEMENT_TYPES) next
        slot_data <- switch(et, points = points, shapes = shapes, NULL)
        if (is.null(slot_data) || !reg %in% names(slot_data)) {
            orphan_inst[[i]] <- TRUE
            next
        }
        layer <- slot_data[[reg]]
        if (is.null(layer) || nrow(layer) == 0L) {
            orphan_inst[[i]] <- TRUE
            next
        }
        layer_ids <- if ("instance_id" %in% colnames(layer))
            as.character(layer[["instance_id"]]) else as.character(seq_len(nrow(layer)))
        if (!as.character(inst) %in% layer_ids)
            orphan_inst[[i]] <- TRUE
    }
    if (any(orphan_inst)) {
        drops[["spatialMap_instance"]] <- spatialMap[orphan_inst, , drop = FALSE]
        spatialMap <- spatialMap[!orphan_inst, , drop = FALSE]
    }
    if (nrow(spatialMap) < n_before && !length(drops))
        drops[["spatialMap"]] <- "rows removed during harmonization"
    list(spatialMap = spatialMap, drops = drops)
}

#' Build a spatialMap table from sampleMap
#'
#' @description
#' Assemble a \code{spatialMap} \linkS4class{DataFrame} linking assay observations
#' to spatial layer instances, analogous to \code{\link[MultiAssayExperiment]{listToMap}}
#' for \code{sampleMap}.
#'
#' @param sampleMap A \linkS4class{DataFrame} with MAE columns \code{assay},
#'   \code{primary}, and \code{colname}.
#' @param region Character. Name of the points or shapes layer (must match a name in
#'   \code{spatialPoints} or \code{spatialShapes} when the object is constructed).
#' @param element_type Character. \code{"points"} or \code{"shapes"}.
#' @param assays Optional character vector limiting which assays from \code{sampleMap}
#'   are included.
#' @param instance_from Character. Column of \code{sampleMap} used for
#'   \code{instance_id} values (default \code{"colname"}, i.e. observation id equals
#'   spatial instance id).
#'
#' @return A \linkS4class{DataFrame} with columns \code{assay}, \code{colname},
#'   \code{element_type}, \code{region}, and \code{instance_id}.
#'
#' @seealso \code{\link{prepMASE}}, \code{\link[MultiAssayExperiment]{listToMap}}
#'
#' @examples
#' sm <- S4Vectors::DataFrame(
#'     assay = factor(rep("rna", 3), "rna"),
#'     primary = c("S1", "S1", "S2"),
#'     colname = paste0("spot", 1:3))
#' buildSpatialMap(sm, region = "tissue1", element_type = "points")
#'
#' @export
buildSpatialMap <-
function(sampleMap, region, element_type = c("points", "shapes"), assays = NULL,
         instance_from = "colname")
{
    element_type <- match.arg(element_type)
    if (!is(sampleMap, "DataFrame"))
        sampleMap <- DataFrame(sampleMap)
    req <- c("assay", "colname")
    if (!all(req %in% colnames(sampleMap)))
        stop(wmsg("'sampleMap' must contain columns: ",
                  paste(req, collapse = ", ")))
    if (!instance_from %in% colnames(sampleMap))
        stop(wmsg("'instance_from' column '", instance_from,
                  "' not found in sampleMap"))
    sm <- sampleMap
    if (!is.null(assays))
        sm <- sm[sm[["assay"]] %in% assays, , drop = FALSE]
    if (nrow(sm) == 0L)
        return(DataFrame(assay = factor(levels = levels(sm[["assay"]])),
               colname = character(),
               element_type = character(),
               region = character(),
               instance_id = character()))
    if (!is.factor(sm[["assay"]]))
        sm[["assay"]] <- factor(sm[["assay"]])
    DataFrame(assay = sm[["assay"]],
              colname = sm[["colname"]],
              element_type = rep(element_type, nrow(sm)),
              region = rep(region, nrow(sm)),
              instance_id = sm[[instance_from]])
}

#' Check spatialMap consistency
#'
#' @description
#' Run \code{spatialMap} diagnostics without constructing a
#' \linkS4class{MultiAssaySpatialExperiment}. Returns a character vector of
#' problems (empty if valid).
#'
#' @param spatialMap A \code{spatialMap} \linkS4class{DataFrame} or \code{NULL}.
#' @param sampleMap A MAE \code{sampleMap} \linkS4class{DataFrame}.
#' @param points A \linkS4class{PointsLayerList} or named list of point layers.
#' @param shapes A \linkS4class{ShapesLayerList} or named list of shape layers.
#'
#' @return Character vector of error messages (length zero if no problems).
#'
#' @examples
#' sm <- S4Vectors::DataFrame(
#'     assay = factor("rna", "rna"),
#'     primary = "S1",
#'     colname = "spot1")
#' spmap <- buildSpatialMap(sm, "coords", "points")
#' spmap[["instance_id"]] <- "missing"
#' checkSpatialMap(spmap, sm, points = PointsLayerList(
#'     coords = S4Vectors::DataFrame(
#'         x = 1, y = 1, instance_id = "spot1")))
#'
#' @export
#' @importFrom S4Vectors DataFrame wmsg
#' @importFrom MultiAssayExperiment ExperimentList
checkSpatialMap <-
function(spatialMap, sampleMap, points = PointsLayerList(),
         shapes = ShapesLayerList())
{
    if (is.null(spatialMap) || nrow(spatialMap) == 0L)
        return(character())
    if (!is(sampleMap, "DataFrame"))
        sampleMap <- DataFrame(sampleMap)
    pts <- if (is(points, "PointsLayerList")) points else PointsLayerList(points)
    shps <- if (is(shapes, "ShapesLayerList")) shapes else ShapesLayerList(shapes)
    primaries <- unique(as.character(sampleMap[["primary"]]))
    assays <- unique(as.character(sampleMap[["assay"]]))
    exps <- lapply(assays, function(a) {
        cn <- unique(as.character(sampleMap[["colname"]][sampleMap[["assay"]] == a]))
        matrix(0, 1L, length(cn), dimnames = list("g1", cn))
    })
    names(exps) <- assays
    tmp <- tryCatch(
        MultiAssaySpatialExperiment(
            experiments = ExperimentList(exps),
            colData = DataFrame(row.names = primaries),
            sampleMap = sampleMap,
            points = pts,
            shapes = shps,
            spatialMap = spatialMap
        ),
        error = function(e) e
    )
    if (inherits(tmp, "error"))
        return(conditionMessage(tmp))
    msg <- tryCatch(validObject(tmp), error = function(e) conditionMessage(e))
    if (isTRUE(msg)) character() else msg
}
