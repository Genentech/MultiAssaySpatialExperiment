### =========================================================================
### MultiAssaySpatialExperiment - SpatialData integration (deferred)
### -------------------------------------------------------------------------
###
### Loaded on demand via .onLoad packageEvent hook when spatialdataR /
### SpatialData is installed. Coercion to/from SpatialData; methods for
### SpatialData generics (images, labels, points, shapes, tables, layer,
### element). Registered in .register_SpatialData_all().
###
### -------------------------------------------------------------------------

#' @include MultiAssaySpatialExperiment-class.R
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion helpers
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @importFrom S4Vectors DataFrame
.pointFrameToDataFrame <- function(pf) {
    dat <- SpatialData::data(pf)
    df <- if (is(dat, "DataFrame")) dat else DataFrame(dat)
    if ("id" %in% colnames(df) && !"instance_id" %in% colnames(df))
        df[["instance_id"]] <- df[["id"]]
    df
}

#' @importFrom S4Vectors DataFrame
.shapeFrameToDataFrame <- function(sf) {
    dat <- SpatialData::data(sf)
    df <- if (is(dat, "DataFrame")) dat else DataFrame(dat)
    if ("id" %in% colnames(df) && !"instance_id" %in% colnames(df))
        df[["instance_id"]] <- df[["id"]]
    df
}

`%||%` <- function(x, y) if (is.null(x)) y else x

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### as(SpatialData, "MultiAssaySpatialExperiment") / as(MASE, "SpatialData")
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Coercion between SpatialData and MultiAssaySpatialExperiment
#'
#' @description
#' \code{as(sd, "MultiAssaySpatialExperiment")} converts a SpatialData object.
#' Tables become assays; points and shapes are collected into the MASE slots.
#'
#' \code{as(mase, "SpatialData")} converts a MASE to SpatialData. Experiments
#' become tables; points and shapes become PointFrame and ShapeFrame elements.
#'
#' @name coerce-SpatialData
#' @aliases coerce,SpatialData,MultiAssaySpatialExperiment-method
#'   coerce,MultiAssaySpatialExperiment,SpatialData-method
NULL

.build_spatialMap_from_sd <- function(tbls, tnames, pts, shps) {
    sda <- "spatialdata_attrs"
    spmap_rows <- list()
    for (k in seq_along(tbls)) {
        sce <- tbls[[k]]
        md <- tryCatch(
            int_metadata(sce)[[sda]],
            error = function(e) NULL
        )
        if (is.null(md) || !is.list(md)) next
        region <- md[["region"]]
        if (is.null(region) || !is.character(region) || length(region) != 1L)
            next
        rk <- md[["region_key"]] %||% "rk"
        ik <- md[["instance_key"]] %||% "ik"
        ass <- tnames[[k]]
        colnms <- colnames(sce)
        if (is.null(colnms) || length(colnms) != ncol(sce))
            colnms <- paste0("cell", seq_len(ncol(sce)))
        if (ik %in% names(colData(sce))) {
            inst <- colData(sce)[[ik]]
        } else {
            inst <- seq_len(ncol(sce))
        }
        
        ## Determine element_type by searching which spatialdata slot contains the region
        ## This aligns with spatialdata's implicit search behavior
        in_pts <- region %in% names(pts)
        in_shps <- region %in% names(shps)
        
        if (in_pts || in_shps) {
            ## Set element_type based on which slot contains the region
            ## If in both (ambiguous), points takes precedence (matching spatialdata behavior)
            element_type <- if (in_pts) "points" else "shapes"
            
            spmap_rows[[length(spmap_rows) + 1L]] <- DataFrame(
                assay = factor(rep(ass, length(colnms)), tnames),
                colname = colnms,
                element_type = rep(element_type, length(colnms)),
                region = rep(region, length(colnms)),
                instance_id = inst
            )
        }
    }
    if (length(spmap_rows) == 0L) return(NULL)
    do.call(rbind, spmap_rows)
}

#' @importClassesFrom MultiAssaySpatialExperiment MultiAssaySpatialExperiment
#' @importFrom MultiAssayExperiment ExperimentList MultiAssayExperiment
#' @importFrom S4Vectors DataFrame
.coerce_SpatialData_to_MASE <- function(from) {
    tbls <- SpatialData::tables(from)
    if (length(tbls) == 0L)
        stop("Coercion from SpatialData requires at least one table; ",
             "found 0")
    tnames <- SpatialData::tableNames(from) %||% paste0("table", seq_along(tbls))
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
    pts <- PointsLayerList()
    pt_names <- SpatialData::pointNames(from)
    if (length(pt_names) > 0L) {
        pt_list <- lapply(pt_names, function(nm) {
            pf <- SpatialData::point(from, nm)
            .pointFrameToDataFrame(pf)
        })
        names(pt_list) <- pt_names
        pts <- PointsLayerList(pt_list)
    }
    shps <- ShapesLayerList()
    shp_names <- SpatialData::shapeNames(from)
    if (length(shp_names) > 0L) {
        shp_list <- lapply(shp_names, function(nm) {
            sf <- SpatialData::shape(from, nm)
            .shapeFrameToDataFrame(sf)
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

#' @importFrom SingleCellExperiment int_metadata int_metadata<-
#' @importFrom SummarizedExperiment colData colData<-
.augment_sce_with_spatialdata_attrs <- function(sce, assay_name, spmap, rk, ik) {
    if (is.null(spmap) || nrow(spmap) == 0L) return(sce)
    sub <- spmap[as.character(spmap[["assay"]]) == assay_name, , drop = FALSE]
    if (nrow(sub) == 0L) return(sce)
    region <- sub[["region"]][[1L]]
    midx <- match(colnames(sce), sub[["colname"]])
    inst_vals <- if (anyNA(midx)) seq_len(ncol(sce)) else sub[["instance_id"]][midx]
    cd <- colData(sce)
    if (!rk %in% names(cd)) cd[[rk]] <- rep(region, nrow(cd))
    if (!ik %in% names(cd)) cd[[ik]] <- inst_vals
    colData(sce) <- cd
    imd <- int_metadata(sce)
    imd[["spatialdata_attrs"]] <- list(region = region, region_key = rk, instance_key = ik)
    int_metadata(sce) <- imd
    sce
}

#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom MultiAssayExperiment experiments
.coerce_MASE_to_SpatialData <- function(from) {
    exps <- experiments(from)
    exp_names <- names(exps) %||% paste0("assay", seq_along(exps))
    spmap <- spatialMap(from)
    rk <- "rk"
    ik <- "instance_id"
    
    ## Note: spatialMap element_type is used during MASE construction/validation
    ## to route regions to correct MASE slots (points vs shapes). When converting
    ## to SpatialData, we simply place each region in its corresponding spatialdata
    ## slot, and the spatialdata_attrs metadata in each table stores only the
    ## region name (not element_type), because spatialdata searches across all
    ## slots to resolve region references.
    
    tbl_list <- lapply(seq_along(exps), function(k) {
        exp <- exps[[k]]
        if (!is(exp, "SingleCellExperiment"))
            exp <- as(exp, "SingleCellExperiment")
        .augment_sce_with_spatialdata_attrs(exp, exp_names[[k]], spmap, rk, ik)
    })
    names(tbl_list) <- exp_names
    pts <- spatialPoints(from)
    pt_names <- spatialPointNames(from)
    pt_list <- if (length(pts) > 0L && length(pt_names) > 0L) {
        lapply(pt_names, function(nm) {
            df <- pts[[nm]]
            SpatialData::PointFrame(data = as.data.frame(df))
        })
    } else {
        list()
    }
    if (length(pt_list) > 0L)
        names(pt_list) <- pt_names
    shps <- spatialShapes(from)
    shp_names <- spatialShapeNames(from)
    shp_list <- if (length(shps) > 0L && length(shp_names) > 0L) {
        lapply(shp_names, function(nm) {
            df <- shps[[nm]]
            SpatialData::ShapeFrame(data = as.data.frame(df))
        })
    } else {
        list()
    }
    if (length(shp_list) > 0L)
        names(shp_list) <- shp_names
    SpatialData::SpatialData(images = list(), labels = list(), points = pt_list,
                             shapes = shp_list, tables = tbl_list)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SpatialData method registration
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.LAYERS <- c("images", "labels", "points", "shapes", "tables")

.register_SpatialData_coercion <- function() {
    setAs("SpatialData", "MultiAssaySpatialExperiment", .coerce_SpatialData_to_MASE)
    setAs("MultiAssaySpatialExperiment", "SpatialData", .coerce_MASE_to_SpatialData)
}

.register_SpatialData_methods <- function() {
    setMethod("images", "MultiAssaySpatialExperiment", function(x) spatialImages(x))
    setMethod("labels", "MultiAssaySpatialExperiment", function(x) spatialLabels(x))
    setMethod("points", "MultiAssaySpatialExperiment", function(x) spatialPoints(x))
    setMethod("shapes", "MultiAssaySpatialExperiment", function(x) spatialShapes(x))
    setMethod("tables", "MultiAssaySpatialExperiment", function(x) experiments(x))

    setMethod("imageNames", "MultiAssaySpatialExperiment", function(x) names(spatialImages(x)))
    setMethod("labelNames", "MultiAssaySpatialExperiment", function(x) names(spatialLabels(x)))
    setMethod("pointNames", "MultiAssaySpatialExperiment", function(x) names(spatialPoints(x)))
    setMethod("shapeNames", "MultiAssaySpatialExperiment", function(x) names(spatialShapes(x)))
    setMethod("tableNames", "MultiAssaySpatialExperiment", function(x) names(experiments(x)))

    setMethod("image", "MultiAssaySpatialExperiment", function(x, i = 1) spatialImages(x)[[i]])
    setMethod("label", "MultiAssaySpatialExperiment", function(x, i = 1) spatialLabels(x)[[i]])
    setMethod("point", "MultiAssaySpatialExperiment", function(x, i = 1) spatialPoints(x)[[i]])
    setMethod("shape", "MultiAssaySpatialExperiment", function(x, i = 1) spatialShapes(x)[[i]])
    ## table conflicts with base::table; target SpatialData's generic explicitly
    tdef <- getGeneric("table", where = asNamespace("SpatialData"))
    if (!is.null(tdef))
        setMethod(tdef, "MultiAssaySpatialExperiment", function(x, i = 1) experiments(x)[[i]])

    setReplaceMethod("images", c("MultiAssaySpatialExperiment", "list"),
        function(x, value) {
            spatialImages(x) <- as(value, "RasterLayerList")
            x
        })
    setReplaceMethod("labels", c("MultiAssaySpatialExperiment", "list"),
        function(x, value) {
            spatialLabels(x) <- as(value, "RasterLayerList")
            x
        })
    setReplaceMethod("points", c("MultiAssaySpatialExperiment", "list"),
        function(x, value) {
            spatialPoints(x) <- as(value, "PointsLayerList")
            x
        })
    setReplaceMethod("shapes", c("MultiAssaySpatialExperiment", "list"),
        function(x, value) {
            spatialShapes(x) <- as(value, "ShapesLayerList")
            x
        })
    setReplaceMethod("tables", c("MultiAssaySpatialExperiment", "list"),
        function(x, value) {
            experiments(x) <- as(value, "ExperimentList")
            x
        })
    setReplaceMethod("tables", c("MultiAssaySpatialExperiment", "List"),
        function(x, value) {
            experiments(x) <- as(value, "ExperimentList")
            x
        })

    .sd_set_one <- function(slot_getter, slot_setter, x, i, value, prefix) {
        y <- slot_getter(x)
        n <- length(y)
        nms <- names(y)
        if (missing(i)) i <- 1L
        if (is.numeric(i) && length(i) == 1L) {
            i <- if (i > n) paste0(prefix, n + 1L) else nms[i]
        }
        if (is.null(value)) {
            y <- y[setdiff(nms, i)]
        } else {
            y[[i]] <- value
        }
        x <- slot_setter(x, y)
        x
    }

    setReplaceMethod("image", c("MultiAssaySpatialExperiment", "ANY", "ANY"),
        function(x, i, value) {
            .sd_set_one(spatialImages, `spatialImages<-`, x, i, value, "image")
        })
    setReplaceMethod("label", c("MultiAssaySpatialExperiment", "ANY", "ANY"),
        function(x, i, value) {
            .sd_set_one(spatialLabels, `spatialLabels<-`, x, i, value, "label")
        })
    setReplaceMethod("point", c("MultiAssaySpatialExperiment", "ANY", "ANY"),
        function(x, i, value) {
            .sd_set_one(spatialPoints, `spatialPoints<-`, x, i, value, "point")
        })
    setReplaceMethod("shape", c("MultiAssaySpatialExperiment", "ANY", "ANY"),
        function(x, i, value) {
            .sd_set_one(spatialShapes, `spatialShapes<-`, x, i, value, "shape")
        })
    setReplaceMethod("table", c("MultiAssaySpatialExperiment", "ANY", "ANY"),
        function(x, i, value) {
            .sd_set_one(experiments, `experiments<-`, x, i, value, "table")
        })

    setMethod("layer", c("MultiAssaySpatialExperiment", "character"),
        function(x, i) {
            i <- match.arg(i, .LAYERS, TRUE)
            switch(i,
                images = spatialImages(x),
                labels = spatialLabels(x),
                points = spatialPoints(x),
                shapes = spatialShapes(x),
                tables = experiments(x))
        })

    setMethod("layer", c("MultiAssaySpatialExperiment", "numeric"),
        function(x, i) {
            if (length(i) != 1L || i <= 0L || i > 5L || i != round(i))
                stop("invalid 'i'; should be an integer in [1, 5], or a ",
                     "string in ", dQuote(paste(.LAYERS, collapse = "/")))
            SpatialData::layer(x, .LAYERS[i])
        })

    setMethod("element", c("MultiAssaySpatialExperiment", "ANY", "character"),
        function(x, i, j) {
            y <- SpatialData::layer(x, i)
            j <- match.arg(j, names(y))
            y[[j]]
        })

    setMethod("element", c("MultiAssaySpatialExperiment", "ANY", "numeric"),
        function(x, i, j) {
            y <- SpatialData::layer(x, i)
            n <- length(y)
            if (n == 0L) stop("there aren't any ", dQuote(i))
            if (is.infinite(j)) j <- n
            if (length(j) != 1L || j <= 0L || j > n || j != round(j))
                stop("invalid 'j'; should be a scalar integer or ",
                     "a string specifying an element in layer 'i'")
            j <- names(y)[j]
            SpatialData::element(x, i, j)
        })
}

.register_SpatialData_activate <- function() {
    .register_SpatialData_coercion()
    .register_SpatialData_methods()
    invisible(NULL)
}
