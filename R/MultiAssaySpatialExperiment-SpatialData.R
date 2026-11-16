### =========================================================================
### MultiAssaySpatialExperiment - spatialdataR integration
### -------------------------------------------------------------------------
###
### Coercion to and from spatialdataR's SpatialData class, and methods for its
### generics (images, labels, points, shapes, tables, layer, element). Nothing
### here is registered until spatialdataR is loaded: .register_SpatialData_all()
### (in zzz.R) imports the generics into .sd_ext_env and registers every method
### there, because the package namespace is sealed by then, as for
### SpatialFeatureExperiment. See ?"coerce-SpatialData".
###
### -------------------------------------------------------------------------

#' @include MultiAssaySpatialExperiment-class.R
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion helpers
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MASE links a table column to a layer row through the layer's "instance_id".
# The column a SpatialData table links on is the table's own instance_key,
# which can have any name, so that column wins; "id" is the fallback.
.sd_with_instance_id <- function(df, key = NULL) {
    if (!is.null(key) && !identical(key, "instance_id") &&
        key %in% colnames(df))
        df[["instance_id"]] <- df[[key]]
    else if ("id" %in% colnames(df) && !"instance_id" %in% colnames(df))
        df[["instance_id"]] <- df[["id"]]
    df
}

# region name -> the instance_key of the (first) table that annotates it
.sd_region_keys <- function(tbls) {
    keys <- character(0L)
    for (k in seq_along(tbls)) {
        md <- int_metadata(tbls[[k]])[["spatialdata_attrs"]]
        if (is.null(md) || !is.list(md)) next
        ik <- md[["instance_key"]]
        if (!is.character(ik) || length(ik) != 1L) next
        for (r in setdiff(as.character(md[["region"]]), names(keys)))
            keys[[r]] <- ik
    }
    keys
}

#' @importFrom S4Vectors DataFrame
.pointFrameToDataFrame <- function(pf, key = NULL) {
    .sd_with_instance_id(DataFrame(as.data.frame(pf), check.names = FALSE),
                         key)
}

#' @importFrom S4Vectors DataFrame
.shapeFrameToDataFrame <- function(sf, key = NULL) {
    .sd_with_instance_id(DataFrame(as.data.frame(sf), check.names = FALSE),
                         key)
}

# The sf that spatialdataR holds for one points or shapes layer. A geometry
# (sfc) column is used as is. Otherwise points take their coordinates from
# "x" and "y" or, as elsewhere in MASE, from the first two columns that are not
# "instance_id".
.sd_frame_data <- function(df, kind, name) {
    d <- as.data.frame(df, optional = TRUE)
    is_sfc <- vapply(d, inherits, NA, what = "sfc")
    if (any(is_sfc))
        return(sf::st_as_sf(d, sf_column_name = names(d)[which(is_sfc)[1L]]))
    if (nrow(d) == 0L) {
        empty <- if (kind == "point") sf::st_point() else sf::st_polygon()
        return(sf::st_sf(geometry = sf::st_sfc(empty)[0L], crs = NA))
    }
    if (kind != "point")
        stop("shapes layer '", name, "' has no geometry column", call. = FALSE)
    cc <- if (all(c("x", "y") %in% names(d))) c("x", "y")
          else setdiff(names(d), "instance_id")[seq_len(2L)]
    if (anyNA(cc) || !all(vapply(d[cc], is.numeric, NA)))
        stop("points layer '", name, "' needs numeric coordinate columns: ",
             "\"x\" and \"y\", or its first two columns other than ",
             "\"instance_id\"", call. = FALSE)
    sf::st_as_sf(d, coords = cc, crs = NA)
}

# The DuckDB connection spatialdataR builds and reads every element on.
.sd_conn <- function() {
    conn <- tryCatch(get(".conn", envir = asNamespace("spatialdataR"))(),
                     error = function(e) NULL)
    if (is.null(conn)) duckspatial::ddbs_default_conn() else conn
}

# spatialdataR's own constructors write every points frame to one fixed table
# name and every shapes frame to another, replacing the table each time, so a
# second element silently overwrites the data of the first. Each element is
# written to a table of its own here instead, and handed over ready-made.
.sd_element_frame <- function(d, kind, tbl) {
    force(d)
    conn <- .sd_conn()
    ## ddbs_write_table() announces the missing CRS and the import; nothing
    ## else is muted
    withCallingHandlers(
        duckspatial::ddbs_write_table(conn = conn, data = d, name = tbl,
                                      overwrite = FALSE, temp_view = FALSE),
        message = function(m) {
            if (grepl("No CRS found|successfully imported",
                      conditionMessage(m)))
                invokeRestart("muffleMessage")
        })
    crs <- sf::st_crs(d)
    duckspatial::as_duckspatial_df(
        x = tbl, conn = conn,
        crs = if (is.na(crs)) NA_character_ else crs$input,
        geom_col = attr(d, "sf_column"))
}

# Drop the tables an interrupted conversion wrote.
.sd_drop_tables <- function(tbls) {
    conn <- .sd_conn()
    for (tbl in tbls)
        try(DBI::dbRemoveTable(conn, tbl), silent = TRUE)
    invisible(NULL)
}

# All elements of one kind, built from a MASE layer list.
.sd_elements <- function(layers, nms, kind, make, written) {
    lapply(stats::setNames(nms, nms), function(nm) {
        ## each step is evaluated here, not passed on as a lazy argument, so an
        ## error in one is raised here, with its own message
        d <- .sd_frame_data(layers[[nm]], kind, nm)
        tbl <- basename(tempfile(paste0("mase_sd", kind, "_")))
        ## recorded before the write, so that a failed write is cleaned up too
        written$tables <- c(written$tables, tbl)
        frame <- .sd_element_frame(d, kind, tbl)
        ## spatialdataR matches a table's instance ids to the element's rows
        ## through the element's instance_key; without one it counts 1..n, and
        ## keeps no row of an element annotated by a table with other ids
        ik <- if ("instance_id" %in% names(d)) "instance_id"
        make(data = frame, ik = ik)
    })
}

# Every element must read back with the rows it was written with. 'get' reads
# an element of the finished SpatialData; spatialdataR keeps only the rows of an
# element that a table annotating it lists, so rows can go missing there.
.sd_check_rows <- function(elements, layers, kind, get = NULL) {
    for (nm in names(elements)) {
        want <- nrow(layers[[nm]])
        got <- length(if (is.null(get)) elements[[nm]] else get(nm))
        if (identical(as.integer(got), as.integer(want))) next
        if (is.null(get))
            stop(kind, " '", nm, "' was written with ", want, " rows but ",
                 "reads back with ", got, "; spatialdataR replaced its data",
                 call. = FALSE)
        warning(kind, " '", nm, "' has ", got, " of its ", want, " rows in ",
                "the SpatialData: spatialdataR keeps only the rows that a ",
                "table annotating it lists by instance id", call. = FALSE)
    }
    invisible(elements)
}

# Points and shapes whose coordinate transformations are anything but identity.
# The coordinates are carried as stored, so those transformations are lost.
.sd_transformed <- function(sd) {
    nms <- c(spatialdataR::pointNames(sd), spatialdataR::shapeNames(sd))
    Filter(function(nm) {
        any(spatialdataR::CTtype(spatialdataR::element(sd, nm)) != "identity")
    }, nms)
}

# A generic of spatialdataR, found in its namespace so that it does not matter
# whether the package is attached. For a base generic that spatialdataR extends,
# such as labels(), this is the generic it derived.
.sd_generic <- function(f) {
    methods::getGeneric(f, where = asNamespace("spatialdataR"))
}

# The generics that MASE adds methods to.
.SD_GENERICS <- c(
    "element",
    "element<-",
    "image",
    "image<-",
    "imageNames",
    "imageNames<-",
    "images",
    "images<-",
    "label",
    "label<-",
    "labelNames",
    "labelNames<-",
    "labels",
    "labels<-",
    "layer",
    "point",
    "point<-",
    "pointNames",
    "pointNames<-",
    "points",
    "points<-",
    "shape",
    "shape<-",
    "shapeNames",
    "shapeNames<-",
    "shapes",
    "shapes<-",
    "table",
    "table<-",
    "tableNames",
    "tableNames<-",
    "tables",
    "tables<-")

# Every method is registered in .sd_ext_env, so the generics are put there too.
.import_SD_generics <- function() {
    gens <- lapply(.SD_GENERICS, .sd_generic)
    lacking <- .SD_GENERICS[vapply(gens, is.null, NA)]
    if (length(lacking))
        stop("this version of spatialdataR lacks the generics ",
             toString(lacking), call. = FALSE)
    for (i in seq_along(gens))
        assign(.SD_GENERICS[[i]], gens[[i]], envir = .sd_ext_env)
    invisible(NULL)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### as(SpatialData, "MultiAssaySpatialExperiment") / as(MASE, "SpatialData")
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Interoperating with spatialdataR
#'
#' @description
#' When the \pkg{spatialdataR} package is loaded, a
#' \code{MultiAssaySpatialExperiment} converts to and from its
#' \code{SpatialData} class and answers its generics. Nothing is registered
#' until \pkg{spatialdataR} is loaded, and it does not matter which package is
#' loaded first.
#'
#' @details
#' \code{as(sd, "MultiAssaySpatialExperiment")} turns each table of the
#' \code{SpatialData} \code{sd} into an assay and collects its points and shapes
#' into the points and shapes layers. A table whose \code{spatialdata_attrs}
#' region names a points or shapes element is linked to that element through
#' \code{spatialMap}. The \code{SpatialData} must hold at least one table.
#'
#' \code{as(mase, "SpatialData")} writes each assay as a table, first coercing
#' it to a \code{SingleCellExperiment} if it is not one, and the points and
#' shapes layers as points and shapes.
#'
#' Points and shapes are loaded into memory, because \pkg{spatialdataR} holds
#' them as lazy DuckDB queries. Going the other way, each element is written to
#' a table of its own on \pkg{spatialdataR}'s DuckDB connection, and its rows
#' are checked once they are all built. Points take their coordinates from the
#' columns \code{x} and \code{y}, or else from the first two columns other than
#' \code{instance_id}, as elsewhere in this package.
#'
#' A table links to a layer row through the table's own \code{instance_key},
#' whatever that column is called in the points or shapes element. The layer
#' gets an \code{instance_id} column copied from it (or from \code{id} when no
#' table names the element).
#'
#' A table's region and instance keys are kept where \pkg{spatialdataR} keeps
#' them (in the table's \code{int_colData}) and read from there, and the CRS of
#' points and shapes is carried over. \pkg{spatialdataR} keeps only the rows of
#' a shapes element that a table annotating it lists; if that drops rows of a
#' layer, a warning says how many.
#'
#' Going to \code{SpatialData}, the sample-level \code{colData} and the sample
#' names have no place in a \code{SpatialData} and are dropped, with a warning
#' if there are \code{colData} columns. Each points or shapes element is written
#' to a table of its own on \pkg{spatialdataR}'s DuckDB connection. A conversion
#' that fails removes the tables it wrote; after a successful one they stay on
#' the connection for the life of the session, because the \code{SpatialData}
#' reads from them.
#'
#' Element names must be unique across points, shapes and tables, because
#' \pkg{spatialdataR} addresses an element by name alone.
#'
#' Images and labels are not converted yet: both directions warn and name what
#' was dropped. So are coordinate transformations: a point or shape element
#' whose transformations are not the identity gets a warning, and its
#' coordinates are used as stored. When converting to \code{SpatialData}, a
#' table whose region annotation names an element that was not converted loses
#' that annotation, with a warning, so that the result is still a valid
#' \code{SpatialData}. A table keeps the region and instance keys it came with.
#' If the \code{spatialMap} covers only some of its columns, the others get an
#' \code{NA} region and instance id, with a warning.
#'
#' @return
#' The coercions return the converted object. The generics return what they do
#' for a \code{SpatialData}, except that \code{element()}, \code{image()},
#' \code{point()} and \code{shape()} return the layer's own object, such as a
#' \code{DataFrame} for a points layer, not a \pkg{spatialdataR} element.
#'
#' @section Generics:
#' With \pkg{spatialdataR} loaded, these of its generics have methods for a
#' \code{MultiAssaySpatialExperiment}: \code{images}, \code{labels},
#' \code{points}, \code{shapes} and \code{tables} (and their replacement forms);
#' \code{imageNames}, \code{labelNames}, \code{pointNames},
#' \code{shapeNames} and
#' \code{tableNames}; and \code{image}, \code{label}, \code{point}, \code{shape}
#' and \code{table} (and their replacement forms).
#'
#' \code{layer(x, i)} returns the layer that holds the element named \code{i}:
#' \code{"images"}, \code{"labels"}, \code{"points"}, \code{"shapes"} or
#' \code{"tables"}. \code{element(x, i)} returns the element itself, given its
#' name or its position counted across the layers in that order, and
#' \code{element(x, i) <- value} replaces it. If a name occurs in more than one
#' layer, \code{element()} stops and says so.
#'
#' Looking up a name or position that is not there is an error for every
#' accessor, as in \pkg{spatialdataR}; an empty layer gives \code{NULL}.
#' Setting an element by position needs one that exists, and a new element needs
#' a name. Assigning \code{NULL} to a layer empties it. The names of the images,
#' labels, points and shapes can be replaced (\code{pointNames(x) <- ...}, and
#' the \code{spatialMap} follows), but not those of the tables, which key the
#' \code{sampleMap} and the \code{spatialMap}; and a table cannot be added
#' through \code{table(x, i) <- value}, which needs a matching \code{sampleMap}.
#'
#' @name coerce-SpatialData
#' @aliases images,MultiAssaySpatialExperiment-method
#' @aliases labels,MultiAssaySpatialExperiment-method
#' @aliases points,MultiAssaySpatialExperiment-method
#' @aliases shapes,MultiAssaySpatialExperiment-method
#' @aliases tables,MultiAssaySpatialExperiment-method
#' @aliases imageNames,MultiAssaySpatialExperiment-method
#' @aliases labelNames,MultiAssaySpatialExperiment-method
#' @aliases pointNames,MultiAssaySpatialExperiment-method
#' @aliases shapeNames,MultiAssaySpatialExperiment-method
#' @aliases tableNames,MultiAssaySpatialExperiment-method
#' @aliases image,MultiAssaySpatialExperiment-method
#' @aliases label,MultiAssaySpatialExperiment-method
#' @aliases point,MultiAssaySpatialExperiment-method
#' @aliases shape,MultiAssaySpatialExperiment-method
#' @aliases table,MultiAssaySpatialExperiment-method
#' @aliases layer,MultiAssaySpatialExperiment-method
#' @aliases element,MultiAssaySpatialExperiment-method
#' @aliases images<-,MultiAssaySpatialExperiment,ANY-method
#' @aliases labels<-,MultiAssaySpatialExperiment,ANY-method
#' @aliases points<-,MultiAssaySpatialExperiment,ANY-method
#' @aliases shapes<-,MultiAssaySpatialExperiment,ANY-method
#' @aliases tables<-,MultiAssaySpatialExperiment,ANY-method
#' @aliases imageNames<-,MultiAssaySpatialExperiment,character-method
#' @aliases labelNames<-,MultiAssaySpatialExperiment,character-method
#' @aliases pointNames<-,MultiAssaySpatialExperiment,character-method
#' @aliases shapeNames<-,MultiAssaySpatialExperiment,character-method
#' @aliases tableNames<-,MultiAssaySpatialExperiment,character-method
#' @aliases image<-,MultiAssaySpatialExperiment,ANY,ANY-method
#' @aliases label<-,MultiAssaySpatialExperiment,ANY,ANY-method
#' @aliases point<-,MultiAssaySpatialExperiment,ANY,ANY-method
#' @aliases shape<-,MultiAssaySpatialExperiment,ANY,ANY-method
#' @aliases table<-,MultiAssaySpatialExperiment,ANY,ANY-method
#' @aliases element<-,MultiAssaySpatialExperiment,character-method
#' @aliases element<-,MultiAssaySpatialExperiment,numeric-method
#' @aliases coerce,SpatialData,MultiAssaySpatialExperiment-method
#' @aliases coerce,MultiAssaySpatialExperiment,SpatialData-method
#'
#' @seealso \code{\link{MultiAssaySpatialExperiment-class}} for the terminology
#' shared with \pkg{spatialdataR}.
#'
#' @examples
#' if (requireNamespace("spatialdataR", quietly = TRUE)) {
#'     zarr <- system.file("extdata", "blobs.zarr", package = "spatialdataR")
#'     sd <- spatialdataR::readSpatialData(zarr)
#'
#'     # images and labels are dropped, with a warning
#'     mase <- suppressWarnings(as(sd, "MultiAssaySpatialExperiment"))
#'     spatialPointNames(mase)
#'     spatialShapeNames(mase)
#'
#'     # the layer that holds an element
#'     spatialdataR::layer(mase, "blobs_points")
#'
#'     sd2 <- suppressWarnings(as(mase, "SpatialData"))
#' }
NULL

#' @importFrom SingleCellExperiment int_colData
.build_spatialMap_from_sd <- function(tbls, tnames, pts, shps) {
    sda <- "spatialdata_attrs"
    spmap_rows <- list()
    for (k in seq_along(tbls)) {
        sce <- tbls[[k]]
        md <- int_metadata(sce)[[sda]]
        if (is.null(md) || !is.list(md)) next
        regions <- as.character(md[["region"]])
        if (!length(regions)) next
        rk <- md[["region_key"]] %||% "rk"
        ik <- md[["instance_key"]] %||% "ik"
        ass <- tnames[[k]]
        colnms <- colnames(sce)
        if (is.null(colnms) || length(colnms) != ncol(sce))
            colnms <- paste0("cell", seq_len(ncol(sce)))
        ## spatialdataR keeps the keys in int_colData; a table built by hand may
        ## have them in colData
        cd <- colData(sce)
        icd <- int_colData(sce)
        key_col <- function(k)
            if (k %in% names(cd)) cd[[k]] else icd[[k]]
        inst <- key_col(ik) %||% seq_len(ncol(sce))

        ## The region each column annotates: its own when the table records it
        ## (a table may annotate several regions), else the single region.
        if (!is.null(key_col(rk))) {
            region <- as.character(key_col(rk))
        } else if (length(regions) == 1L) {
            region <- rep(regions, ncol(sce))
        } else {
            next
        }

        ## The element type comes from the slot that holds the region; points
        ## win when a name is in both, as in spatialdata.
        in_pts <- region %in% names(pts)
        in_shps <- region %in% names(shps)
        keep <- !is.na(region) & !is.na(inst) & (in_pts | in_shps)
        if (!any(keep)) next
        spmap_rows[[length(spmap_rows) + 1L]] <- DataFrame(
            assay = factor(rep(ass, sum(keep)), tnames),
            colname = colnms[keep],
            element_type = ifelse(in_pts[keep], "points", "shapes"),
            region = region[keep],
            instance_id = inst[keep]
        )
    }
    if (length(spmap_rows) == 0L) return(NULL)
    do.call(rbind, spmap_rows)
}

#' @importFrom MultiAssayExperiment ExperimentList MultiAssayExperiment
#' @importFrom S4Vectors DataFrame
.coerce_SpatialData_to_MASE <- function(from) {
    tbls <- spatialdataR::tables(from)
    if (length(tbls) == 0L)
        stop("Coercion from SpatialData requires at least one table; ",
             "found 0")
    dropped <- c(spatialdataR::imageNames(from), spatialdataR::labelNames(from))
    if (length(dropped))
        warning("images and labels are not converted and are dropped: ",
                paste(dropped, collapse = ", "), call. = FALSE)
    moved <- .sd_transformed(from)
    if (length(moved))
        warning("coordinate transformations are not converted; the ",
                "coordinates of these are used as stored: ",
                paste(moved, collapse = ", "), call. = FALSE)
    tnames <- spatialdataR::tableNames(from) %||%
        paste0("table", seq_along(tbls))
    experiments <- ExperimentList(setNames(tbls, tnames))
    primaries <- tnames
    cd <- DataFrame(row.names = primaries)
    sm_rows <- lapply(seq_along(tbls), function(k) {
        sce <- tbls[[k]]
        pname <- tnames[[k]]
        colnms <- colnames(sce)
        if (is.null(colnms) || length(colnms) != ncol(sce))
            colnms <- paste0("cell", seq_len(ncol(sce)))
        n <- ncol(sce)
        DataFrame(assay = factor(rep(pname, n), tnames),
                  primary = rep(pname, n),
                  colname = colnms)
    })
    sm <- do.call(rbind, sm_rows)
    mae <- MultiAssayExperiment(experiments = experiments, colData = cd,
                                sampleMap = sm)
    rkeys <- .sd_region_keys(tbls)
    pts <- PointsLayerList()
    pt_names <- spatialdataR::pointNames(from)
    if (length(pt_names) > 0L) {
        pt_list <- lapply(pt_names, function(nm) {
            pf <- spatialdataR::point(from, nm)
            .pointFrameToDataFrame(pf, unname(rkeys[nm]))
        })
        names(pt_list) <- pt_names
        pts <- PointsLayerList(pt_list)
    }
    shps <- ShapesLayerList()
    shp_names <- spatialdataR::shapeNames(from)
    if (length(shp_names) > 0L) {
        shp_list <- lapply(shp_names, function(nm) {
            sf <- spatialdataR::shape(from, nm)
            .shapeFrameToDataFrame(sf, unname(rkeys[nm]))
        })
        names(shp_list) <- shp_names
        shps <- ShapesLayerList(shp_list)
    }
    spmap <- .build_spatialMap_from_sd(tbls, tnames, pts, shps)
    new("MultiAssaySpatialExperiment",
        mae,
        images = RasterLayerList(),
        labels = RasterLayerList(),
        points = pts,
        shapes = shps,
        imgData = NULL,
        spatialMap = spmap)
}

# A table's spatialdata_attrs$region must name an element that exists in the
# SpatialData, or spatialdataR rejects the object. Images and labels are not
# converted, so a table that annotates one loses that annotation.
#' @importFrom SingleCellExperiment int_metadata int_metadata<-
.sd_drop_stale_region <- function(sce, elements, table) {
    imd <- int_metadata(sce)
    region <- imd[["spatialdata_attrs"]][["region"]]
    if (is.null(region) || all(region %in% elements))
        return(sce)
    warning("table '", table, "' annotates ",
            toString(setdiff(region, elements)), ", which is not converted; ",
            "its region annotation is removed", call. = FALSE)
    imd[["spatialdata_attrs"]] <- NULL
    int_metadata(sce) <- imd
    sce
}

#' @importFrom SingleCellExperiment int_metadata int_metadata<- int_colData
#' @importFrom SingleCellExperiment int_colData<-
#' @importFrom SummarizedExperiment colData colData<-
.augment_sce_with_spatialdata_attrs <-
function(sce, assay_name, spmap, rk, ik)
{
    if (is.null(spmap) || nrow(spmap) == 0L) return(sce)
    sub <- spmap[as.character(spmap[["assay"]]) == assay_name, , drop = FALSE]
    if (nrow(sub) == 0L) return(sce)
    imd <- int_metadata(sce)
    ## a table that came from a SpatialData keeps the keys it had
    old <- imd[["spatialdata_attrs"]]
    rk <- old[["region_key"]] %||% rk
    ik <- old[["instance_key"]] %||% ik
    midx <- match(colnames(sce), sub[["colname"]])
    if (anyDuplicated(sub[["colname"]]))
        warning("table '", assay_name, "': some columns map to more than one ",
                "instance; only the first mapping of each is written",
                call. = FALSE)
    if (anyNA(midx))
        warning("table '", assay_name, "': ", sum(is.na(midx)), " of ",
                length(midx), " columns are not in the spatialMap and get an ",
                "NA region and instance id", call. = FALSE)
    region <- as.character(sub[["region"]])[midx]
    ## spatialdataR keeps the region and instance keys in int_colData, not in
    ## colData, and reads them from there
    icd <- int_colData(sce)
    icd[[rk]] <- region
    icd[[ik]] <- sub[["instance_id"]][midx]
    int_colData(sce) <- icd
    cd <- colData(sce)
    cd[[rk]] <- NULL
    cd[[ik]] <- NULL
    colData(sce) <- cd
    imd[["spatialdata_attrs"]] <- list(
        region = unique(region[!is.na(region)]),
        region_key = rk, instance_key = ik)
    int_metadata(sce) <- imd
    sce
}

#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom MultiAssayExperiment experiments
.coerce_MASE_to_SpatialData <- function(from) {
    exps <- experiments(from)
    exp_names <- names(exps) %||% paste0("assay", seq_along(exps))
    spmap <- spatialMap(from)
    dropped <- c(names(spatialImages(from)), names(spatialLabels(from)))
    if (length(dropped))
        warning("images and labels are not converted and are dropped: ",
                paste(dropped, collapse = ", "), call. = FALSE)
    cd <- colData(from)
    if (ncol(cd) > 0L)
        warning("sample-level colData is not stored in a SpatialData and is ",
                "dropped: ", toString(colnames(cd)), call. = FALSE)
    pts <- spatialPoints(from)
    pt_names <- spatialPointNames(from)
    shps <- spatialShapes(from)
    shp_names <- spatialShapeNames(from)
    all_names <- c(pt_names, shp_names, exp_names)
    if (anyDuplicated(all_names))
        stop("spatialdataR needs element names that are unique across points, ",
             "shapes and tables, but these occur more than once: ",
             toString(unique(all_names[duplicated(all_names)])), call. = FALSE)
    elements <- c(pt_names, shp_names)
    rk <- "rk"
    ik <- "instance_id"

    ## Note: spatialMap element_type is used during MASE construction and
    ## validation to route regions to the right MASE slots (points vs shapes).
    ## Converting to SpatialData, each region is placed in its own spatialdata
    ## slot, and the spatialdata_attrs metadata in each table stores only the
    ## region name (not element_type), because spatialdata searches across all
    ## slots to resolve region references.

    tbl_list <- lapply(seq_along(exps), function(k) {
        exp <- exps[[k]]
        if (!is(exp, "SingleCellExperiment"))
            exp <- as(exp, "SingleCellExperiment")
        exp <- .augment_sce_with_spatialdata_attrs(
            exp, exp_names[[k]], spmap, rk, ik)
        .sd_drop_stale_region(exp, elements, exp_names[[k]])
    })
    names(tbl_list) <- exp_names
    written <- new.env(parent = emptyenv())
    written$tables <- character(0L)
    ok <- FALSE
    on.exit(if (!ok) .sd_drop_tables(written$tables), add = TRUE)
    pt_list <- .sd_elements(pts, pt_names, "point", .sd_make_point, written)
    shp_list <- .sd_elements(shps, shp_names, "shape", .sd_make_shape, written)
    .sd_check_rows(pt_list, pts, "point")
    .sd_check_rows(shp_list, shps, "shape")
    sd <- spatialdataR::SpatialData(images = list(), labels = list(),
                                    points = pt_list, shapes = shp_list,
                                    tables = tbl_list)
    .sd_check_rows(pt_list, pts, "point",
                   function(nm) spatialdataR::point(sd, nm))
    .sd_check_rows(shp_list, shps, "shape",
                   function(nm) spatialdataR::shape(sd, nm))
    ok <- TRUE
    sd
}

.sd_make_point <- function(data, ik = NULL)
    spatialdataR::SpatialDataPoint(data = data, ik = ik)

.sd_make_shape <- function(data, ik = NULL) {
    meta <- spatialdataR::SpatialDataAttrs(type = "shape")
    if (!is.null(ik)) {
        za <- as.list(meta)
        za$spatialdata_attrs <- list(instance_key = ik)
        meta <- spatialdataR::SpatialDataAttrs(za)
    }
    spatialdataR::SpatialDataShape(data = data, meta = meta)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SpatialData method registration
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.LAYERS <- c("images", "labels", "points", "shapes", "tables")

# spatialdataR addresses an element by name, or by position counted across the
# layers; layer(x, i) is the layer that holds element 'i'.
.sd_layer_index <- function(x) {
    list(images = names(spatialImages(x)),
         labels = names(spatialLabels(x)),
         points = names(spatialPoints(x)),
         shapes = names(spatialShapes(x)),
         tables = names(experiments(x)))
}

.sd_layer_getter <- function(l) {
    switch(l, images = spatialImages, labels = spatialLabels,
           points = spatialPoints, shapes = spatialShapes,
           tables = experiments)
}

.sd_layer_setter <- function(l) {
    switch(l, images = `spatialImages<-`, labels = `spatialLabels<-`,
           points = `spatialPoints<-`, shapes = `spatialShapes<-`,
           tables = `experiments<-`)
}

.sd_one_layer <- function(x, i) {
    l <- spatialdataR::layer(x, i)
    if (length(l) > 1L)
        stop("'i' names elements in several layers: ",
             paste(l, collapse = ", "))
    l
}

# One element of a layer, by name or position. An empty layer gives NULL, and
# a name or position that is not there is an error, as in spatialdataR.
.sd_get_one <- function(y, i, what) {
    if (!length(y)) return(NULL)
    if (length(i) != 1L)
        stop("invalid 'i'; should be a single name or position", call. = FALSE)
    if (is.numeric(i)) {
        if (!is.finite(i) || i < 1L || i != round(i))
            stop("invalid 'i'; should be a positive integer or string",
                 call. = FALSE)
        if (i > length(y))
            stop("invalid 'i'; only ", length(y), " ", what,
                 " element(s) available", call. = FALSE)
        i <- names(y)[i]
    }
    if (!i %in% names(y))
        stop("invalid 'i'; should be one of: ", toString(names(y)),
             call. = FALSE)
    y[[i]]
}

# Set, add or (with a NULL value) drop one element of a layer. A position
# must name an element that exists; a new element needs a name.
.sd_set_one <- function(slot_getter, slot_setter, x, i, value, prefix) {
    y <- slot_getter(x)
    nms <- names(y)
    if (missing(i)) i <- 1L
    if (length(i) != 1L)
        stop("invalid 'i'; should be a single name or position", call. = FALSE)
    if (is.numeric(i)) {
        if (!is.finite(i) || i < 1L || i > length(y) || i != round(i))
            stop("invalid 'i'; only ", length(y), " ", prefix,
                 " element(s) available; add a new one by name",
                 call. = FALSE)
        i <- nms[i]
    }
    if (is.null(value)) {
        y <- y[setdiff(nms, i)]
    } else {
        y[[i]] <- value
    }
    slot_setter(x, y)
}

# A whole layer from 'value'; NULL or an empty list empties the layer, as in
# spatialdataR.
.sd_as_layer <- function(value, class) {
    if (is.null(value) || (is.list(value) && !length(value)))
        new(class)
    else
        as(value, class)
}

# Rename the elements of a points, shapes, images or labels layer. The
# spatialMap names regions, so it follows the points and shapes.
.sd_rename_layer <- function(x, layer, value, getter, setter, type = NULL) {
    y <- getter(x)
    if (!is.character(value) || length(value) != length(y) ||
        anyNA(value) || anyDuplicated(value))
        stop("'value' should be ", length(y), " unique names", call. = FALSE)
    old <- names(y)
    others <- setdiff(unlist(.sd_layer_index(x), use.names = FALSE), old)
    if (any(value %in% others))
        stop("spatialdataR needs element names that are unique across layers; ",
             "already used: ", toString(intersect(value, others)),
             call. = FALSE)
    names(y) <- value
    sm <- spatialMap(x)
    x <- setter(x, y)
    if (!is.null(type) && !is.null(sm) && nrow(sm) > 0L) {
        hit <- sm[["element_type"]] == type & sm[["region"]] %in% old
        sm[["region"]][hit] <- value[match(sm[["region"]][hit], old)]
        spatialMap(x) <- sm
    }
    x
}

.register_SpatialData_coercion <- function() {
    setAs("SpatialData", "MultiAssaySpatialExperiment",
          .coerce_SpatialData_to_MASE, where = .sd_ext_env)
    setAs("MultiAssaySpatialExperiment", "SpatialData",
          .coerce_MASE_to_SpatialData, where = .sd_ext_env)
}

.register_SD_accessors <- function() {
    cls <- "MultiAssaySpatialExperiment"
    setMethod("images", cls, function(x) spatialImages(x),
              where = .sd_ext_env)
    ## labels() is base::labels(object, ...), not a spatialdataR-specific
    ## generic
    setMethod("labels", cls, function(object, ...) spatialLabels(object),
              where = .sd_ext_env)
    setMethod("points", cls, function(x) spatialPoints(x),
              where = .sd_ext_env)
    setMethod("shapes", cls, function(x) spatialShapes(x),
              where = .sd_ext_env)
    setMethod("tables", cls, function(x) experiments(x),
              where = .sd_ext_env)

    setMethod("imageNames", cls, function(x) names(spatialImages(x)),
              where = .sd_ext_env)
    setMethod("labelNames", cls, function(x) names(spatialLabels(x)),
              where = .sd_ext_env)
    setMethod("pointNames", cls, function(x) names(spatialPoints(x)),
              where = .sd_ext_env)
    setMethod("shapeNames", cls, function(x) names(spatialShapes(x)),
              where = .sd_ext_env)
    setMethod("tableNames", cls, function(x) names(experiments(x)),
              where = .sd_ext_env)

    setMethod("image", cls, function(x, i = 1) .sd_get_one(spatialImages(x), i, "image"),
              where = .sd_ext_env)
    setMethod("label", cls, function(x, i = 1) .sd_get_one(spatialLabels(x), i, "label"),
              where = .sd_ext_env)
    setMethod("point", cls, function(x, i = 1) .sd_get_one(spatialPoints(x), i, "point"),
              where = .sd_ext_env)
    setMethod("shape", cls, function(x, i = 1) .sd_get_one(spatialShapes(x), i, "shape"),
              where = .sd_ext_env)
    setMethod("table", cls, function(x, i = 1) .sd_get_one(experiments(x), i, "table"),
              where = .sd_ext_env)
}

.register_SD_setters <- function() {
    cls <- "MultiAssaySpatialExperiment"
    setReplaceMethod("images", c(cls, "ANY"), function(x, value) {
        spatialImages(x) <- .sd_as_layer(value, "RasterLayerList")
        x
    }, where = .sd_ext_env)
    setReplaceMethod("labels", c(cls, "ANY"), function(x, value) {
        spatialLabels(x) <- .sd_as_layer(value, "RasterLayerList")
        x
    }, where = .sd_ext_env)
    setReplaceMethod("points", c(cls, "ANY"), function(x, value) {
        spatialPoints(x) <- .sd_as_layer(value, "PointsLayerList")
        x
    }, where = .sd_ext_env)
    setReplaceMethod("shapes", c(cls, "ANY"), function(x, value) {
        spatialShapes(x) <- .sd_as_layer(value, "ShapesLayerList")
        x
    }, where = .sd_ext_env)
    setReplaceMethod("tables", c(cls, "ANY"), function(x, value) {
        experiments(x) <- as(value, "ExperimentList")
        x
    }, where = .sd_ext_env)

    nm_sig <- c(cls, "character")
    setReplaceMethod("imageNames", nm_sig, function(x, value)
        .sd_rename_layer(x, "images", value, spatialImages,
                         `spatialImages<-`), where = .sd_ext_env)
    setReplaceMethod("labelNames", nm_sig, function(x, value)
        .sd_rename_layer(x, "labels", value, spatialLabels,
                         `spatialLabels<-`), where = .sd_ext_env)
    setReplaceMethod("pointNames", nm_sig, function(x, value)
        .sd_rename_layer(x, "points", value, spatialPoints,
                         `spatialPoints<-`, "points"), where = .sd_ext_env)
    setReplaceMethod("shapeNames", nm_sig, function(x, value)
        .sd_rename_layer(x, "shapes", value, spatialShapes,
                         `spatialShapes<-`, "shapes"), where = .sd_ext_env)
    setReplaceMethod("tableNames", nm_sig, function(x, value)
        stop("tables cannot be renamed here: their names key the sampleMap ",
             "and the spatialMap; rebuild the object with the new names",
             call. = FALSE), where = .sd_ext_env)

    sig <- c(cls, "ANY", "ANY")
    setReplaceMethod("image", sig, function(x, i, value) {
        .sd_set_one(spatialImages, `spatialImages<-`, x, i, value, "image")
    }, where = .sd_ext_env)
    setReplaceMethod("label", sig, function(x, i, value) {
        .sd_set_one(spatialLabels, `spatialLabels<-`, x, i, value, "label")
    }, where = .sd_ext_env)
    setReplaceMethod("point", sig, function(x, i, value) {
        .sd_set_one(spatialPoints, `spatialPoints<-`, x, i, value, "point")
    }, where = .sd_ext_env)
    setReplaceMethod("shape", sig, function(x, i, value) {
        .sd_set_one(spatialShapes, `spatialShapes<-`, x, i, value, "shape")
    }, where = .sd_ext_env)
    setReplaceMethod("table", sig, function(x, i, value) {
        .sd_set_one(experiments, `experiments<-`, x, i, value, "table")
    }, where = .sd_ext_env)
}

.register_SD_layer_element <- function() {
    cls <- "MultiAssaySpatialExperiment"
    setMethod("layer", c(cls, "character"), function(x, i, ...) {
        idx <- .sd_layer_index(x)
        if (length(i) != 1L || !i %in% unlist(idx, use.names = FALSE))
            stop("invalid 'i'; should be the name of an element in one of ",
                 paste(.LAYERS, collapse = ", "))
        names(Filter(function(e) i %in% e, idx))
    }, where = .sd_ext_env)
    setMethod("layer", c(cls, "ANY"), function(x, i, ...) {
        stop("invalid 'i'; should be the name of an element")
    }, where = .sd_ext_env)

    setMethod("element", c(cls, "character"), function(x, i, ...) {
        .sd_layer_getter(.sd_one_layer(x, i))(x)[[i]]
    }, where = .sd_ext_env)
    setMethod("element", c(cls, "numeric"), function(x, i, ...) {
        nms <- unlist(.sd_layer_index(x), use.names = FALSE)
        if (length(i) != 1L || i < 1L || i > length(nms) || i != round(i))
            stop("invalid 'i'; should be an integer in [1, ",
                 length(nms), "], or the name of an element")
        spatialdataR::element(x, nms[i])
    }, where = .sd_ext_env)
    setMethod("element", c(cls, "missing"),
              function(x, i, ...) spatialdataR::element(x, 1L),
              where = .sd_ext_env)
    setMethod("element", c(cls, "ANY"), function(x, i, ...) {
        stop("invalid 'i'; should be the name of an element")
    }, where = .sd_ext_env)
    setReplaceMethod("element", c(cls, "character"), function(x, i, value) {
        l <- .sd_one_layer(x, i)
        y <- .sd_layer_getter(l)(x)
        y[[i]] <- value
        .sd_layer_setter(l)(x, y)
    }, where = .sd_ext_env)
    setReplaceMethod("element", c(cls, "numeric"), function(x, i, value) {
        nms <- unlist(.sd_layer_index(x), use.names = FALSE)
        if (length(i) != 1L || i < 1L || i > length(nms) || i != round(i))
            stop("invalid 'i'; should be an integer in [1, ",
                 length(nms), "], or the name of an element", call. = FALSE)
        spatialdataR::element(x, nms[i]) <- value
        x
    }, where = .sd_ext_env)
}

.register_SpatialData_methods <- function() {
    .register_SD_accessors()
    .register_SD_setters()
    .register_SD_layer_element()
}
