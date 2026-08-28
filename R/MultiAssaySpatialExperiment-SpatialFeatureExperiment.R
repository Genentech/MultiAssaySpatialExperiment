### =========================================================================
### MultiAssaySpatialExperiment - SpatialFeatureExperiment integration
### -------------------------------------------------------------------------
###
### Coercion to SpatialFeatureExperiment; methods for dimGeometry, annotGeometry,
### spatialGraphs, bbox, unit, localResults, splitByCol, aggregate, findDebrisCells.
### Registered at runtime via packageEvent when SpatialFeatureExperiment
### (Suggests) is loaded. No NAMESPACE import of SFE classes.
###
### The stored experiment is always a plain SingleCellExperiment (never an SFE
### or SPE) so that all spatial content lives in MASE's own points/shapes/
### spatialMap/imgData layers -- the same layers writeParquet() understands.
### colGeometries beyond the primary centroids become named 'shapes' entries
### with matching spatialMap rows; annotGeometries become 'shapes' entries
### with NO spatialMap rows (that absence is what distinguishes them on the
### way back). Col-margin spatialGraphs are stored as colPair()s on the
### stored SingleCellExperiment, reusing existing SCE/BiocDuckDB machinery
### instead of a new MASE slot. localResults are not carried over (see
### localResults()/localResultNames() below). The 13 SFE-delegating methods
### rebuild an ephemeral SFE on demand via .materialize_SFE_view() (which is
### exactly the MASE -> SpatialFeatureExperiment coercion) rather than reading
### a persisted one.
###
### -------------------------------------------------------------------------

#' @include MultiAssaySpatialExperiment-class.R
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helpers: sf <-> MASE-native DataFrame, listw <-> SelfHits
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @importFrom S4Vectors DataFrame
.sf_to_shapes_df <- function(sfobj, instance_ids = NULL) {
    df <- DataFrame(sfobj)
    if (!is.null(instance_ids))
        df[["instance_id"]] <- instance_ids
    df
}

# Converts a MASE points/shapes DataFrame (a 'geometry' sfc column, or plain
# numeric coordinate columns) back to an sf object. When 'cn' is given and the
# entry has an 'instance_id' column that covers it, rows are reordered to
# match 'cn'.
.shapes_df_to_sf <- function(df, cn = NULL) {
    if (is.null(df) || nrow(df) == 0L)
        return(NULL)
    d <- as.data.frame(df, optional = TRUE)
    if ("instance_id" %in% colnames(d))
        rownames(d) <- as.character(d[["instance_id"]])
    if ("geometry" %in% colnames(d) && inherits(d[["geometry"]], "sfc")) {
        sfobj <- sf::st_as_sf(d)
    } else {
        coord_cols <- setdiff(colnames(d), "instance_id")
        if (length(coord_cols) < 2L)
            return(NULL)
        sfobj <- sf::st_as_sf(d, coords = coord_cols[seq_len(2L)])
    }
    if (!is.null(cn) && !is.null(rownames(sfobj)) && all(cn %in% rownames(sfobj)))
        sfobj <- sfobj[cn, , drop = FALSE]
    sfobj
}

# Lossless for connectivity + weight (see plan doc); 'style' and node names
# are optional/informational and stored in metadata(), not in the edge list
# itself, so they only survive an in-memory round trip, not a Parquet one.
#' @importFrom S4Vectors SelfHits mcols<- metadata<-
.listw_to_selfhits <- function(g) {
    nb <- g[["neighbours"]]
    w <- g[["weights"]]
    nnode <- length(nb)
    from <- rep(seq_len(nnode), lengths(nb))
    to <- unlist(nb, use.names = FALSE)
    weight <- unlist(w, use.names = FALSE)
    hits <- SelfHits(from, to, nnode = nnode)
    mcols(hits) <- DataFrame(weight = weight)
    metadata(hits) <- list(style = g[["style"]])
    hits
}

#' @importFrom S4Vectors from to nnode mcols metadata
.selfhits_to_listw <- function(hits, region_id = NULL) {
    if (!requireNamespace("spdep", quietly = TRUE))
        stop("package 'spdep' is required to materialize a spatial graph ",
             "as a 'listw' object")
    n <- nnode(hits)
    fr <- factor(from(hits), levels = seq_len(n))
    nb <- lapply(split(to(hits), fr), as.integer)
    class(nb) <- "nb"
    if (!is.null(region_id))
        attr(nb, "region.id") <- region_id
    w <- mcols(hits)[["weight"]]
    wl <- unname(split(w, fr))
    style <- metadata(hits)[["style"]]
    if (is.null(style)) style <- "W"
    spdep::nb2listw(nb, glist = wl, style = style, zero.policy = TRUE)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion (registered in .onLoad when SFE available)
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @importClassesFrom SpatialExperiment SpatialExperiment
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom MultiAssayExperiment ExperimentList MultiAssayExperiment experiments
#' @importFrom SpatialExperiment imgData imgData<-
#' @importFrom SingleCellExperiment colPair<-
#' @importFrom S4Vectors DataFrame metadata<-
.register_SFE_coercion <- function() {
    setAs("SpatialFeatureExperiment", "MultiAssaySpatialExperiment",
        function(from) {
            assay_name <- "spatial"
            n <- ncol(from)
            cn <- colnames(from)
            sce <- as(from, "SingleCellExperiment")

            # Primary centroid geometry -> points("coordinates"), exactly
            # like the SpatialExperiment forward coercion.
            scoords <- SpatialExperiment::spatialCoords(from)
            pt_df <- DataFrame(scoords)
            rn <- rownames(scoords)
            pt_df[["instance_id"]] <- if (!is.null(rn) && length(rn) > 0L) rn else cn
            pts <- PointsLayerList(coordinates = pt_df)

            sm <- DataFrame(assay = factor(rep(assay_name, n), assay_name),
                            colname = cn)
            sm_pts <- sm
            sm_pts[[".instance_id"]] <- pt_df[["instance_id"]]
            spatial_map <- buildSpatialMap(sm_pts, region = "coordinates",
                element_type = "points", instance_from = ".instance_id")

            # Every other colGeometry -> a per-cell 'shapes' entry, with
            # matching spatialMap rows.
            shapes_list <- list()
            other_cg <- setdiff(SpatialFeatureExperiment::colGeometryNames(from),
                                "centroids")
            for (nm in other_cg) {
                cg <- SpatialFeatureExperiment::colGeometry(from, type = nm)
                cg_inst <- rownames(cg)
                if (is.null(cg_inst) || length(cg_inst) != nrow(cg))
                    cg_inst <- cn
                shapes_list[[nm]] <- .sf_to_shapes_df(cg, cg_inst)
                sm_cg <- sm
                sm_cg[[".instance_id"]] <- cg_inst
                spatial_map <- .rbindSpatialMaps(spatial_map,
                    buildSpatialMap(sm_cg, region = nm, element_type = "shapes",
                        instance_from = ".instance_id"))
            }

            # annotGeometries -> 'shapes' entries with NO spatialMap rows;
            # that absence is what marks them as annotation-level on the way
            # back (they are not one-row-per-cell).
            ag_names <- SpatialFeatureExperiment::annotGeometryNames(from)
            shapes_list[ag_names] <- lapply(ag_names, function(nm)
                .sf_to_shapes_df(SpatialFeatureExperiment::annotGeometry(from, type = nm)))

            # Col-margin spatial graphs -> colPair()s on the stored SCE.
            graph_names <- tryCatch(SpatialFeatureExperiment::colGraphNames(from),
                error = function(e) character(0L))
            for (nm in graph_names) {
                g <- tryCatch(SpatialFeatureExperiment::colGraph(from, type = nm),
                    error = function(e) NULL)
                if (!is.null(g))
                    colPair(sce, nm) <- .listw_to_selfhits(g)
            }

            lr <- tryCatch(SpatialFeatureExperiment::localResults(from),
                error = function(e) NULL)
            if (!is.null(lr) && length(lr))
                message("'localResults' on the SpatialFeatureExperiment are ",
                        "not carried over by this coercion; access them on ",
                        "the original object before coercing.")

            experiments <- ExperimentList(stats::setNames(list(sce), assay_name))
            sids <- if ("sample_id" %in% colnames(SummarizedExperiment::colData(from)))
                as.character(SummarizedExperiment::colData(from)[["sample_id"]])
            else rep("sample01", n)
            primaries <- unique(sids)
            if (length(primaries) == 0L)
                primaries <- "sample01"
            cd <- DataFrame(row.names = primaries)
            sm_full <- DataFrame(assay = factor(rep(assay_name, n), assay_name),
                                 primary = sids, colname = cn)
            mae <- MultiAssayExperiment(experiments = experiments, colData = cd,
                                        sampleMap = sm_full)

            mase <- new("MultiAssaySpatialExperiment", mae,
                images = RasterLayerList(), labels = RasterLayerList(),
                points = pts, shapes = ShapesLayerList(shapes_list),
                imgData = imgData(from), spatialMap = spatial_map)

            u <- tryCatch(SpatialFeatureExperiment::unit(from),
                error = function(e) NULL)
            if (!is.null(u))
                metadata(mase)[["spatialUnit"]] <- u
            mase
        },
        where = .sfe_ext_env
    )
    setAs("MultiAssaySpatialExperiment", "SpatialFeatureExperiment",
        function(from) .materialize_SFE_view(from),
        where = .sfe_ext_env
    )
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper: materialize an SFE view from MASE-native layers
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @importClassesFrom SpatialExperiment SpatialExperiment
#' @importFrom MultiAssayExperiment experiments
#' @importFrom SingleCellExperiment colPairs
#' @importFrom SummarizedExperiment assays colData rowData
#' @importFrom S4Vectors metadata
.materialize_SFE_view <- function(x, assay = 1L) {
    exps <- experiments(x)
    if (length(exps) != 1L)
        stop("Coercion to SpatialFeatureExperiment requires exactly ",
             "one assay; found ", length(exps))
    exp <- exps[[assay]]
    if (is(exp, "SpatialFeatureExperiment"))
        return(exp)
    if (!requireNamespace("SpatialFeatureExperiment", quietly = TRUE))
        stop("package 'SpatialFeatureExperiment' is required to ",
             "materialize a SpatialFeatureExperiment view")

    cn <- colnames(exp)
    centroids_sf <- .shapes_df_to_sf(spatialPoints(x)[["coordinates"]], cn)
    if (is.null(centroids_sf))
        stop("Cannot materialize a SpatialFeatureExperiment from 'x': ",
             "no resolvable 'coordinates' points layer")

    spmap <- spatialMap(x)
    shp_names <- spatialShapeNames(x)
    per_cell <- character(0L)
    if (!is.null(spmap) && nrow(spmap) > 0L) {
        shp_map <- spmap[spmap[["element_type"]] == "shapes", , drop = FALSE]
        per_cell <- intersect(shp_names, unique(as.character(shp_map[["region"]])))
    }
    annot_names <- setdiff(shp_names, per_cell)

    extra_cgs <- stats::setNames(
        lapply(per_cell, function(nm) .shapes_df_to_sf(spatialShapes(x)[[nm]], cn)),
        per_cell)
    extra_cgs <- Filter(Negate(is.null), extra_cgs)
    annots <- stats::setNames(
        lapply(annot_names, function(nm) .shapes_df_to_sf(spatialShapes(x)[[nm]])),
        annot_names)
    annots <- Filter(Negate(is.null), annots)

    cps <- colPairs(exp)
    graphs <- stats::setNames(
        lapply(names(cps), function(nm) .selfhits_to_listw(cps[[nm]], cn)),
        names(cps))

    u <- metadata(x)[["spatialUnit"]]
    if (is.null(u)) u <- "full_res_image_pixel"

    sample_ids <- unique(as.character(colData(exp)[["sample_id"]]))
    if (length(sample_ids) == 0L)
        sample_ids <- "sample01"
    graphs_arg <- if (length(graphs))
            stats::setNames(lapply(sample_ids, function(sid)
                list(row = NULL, col = graphs, annot = NULL)), sample_ids)
        else NULL

    sfe <- SpatialFeatureExperiment::SpatialFeatureExperiment(
        assays = assays(exp), rowData = rowData(exp), colData = colData(exp),
        colGeometries = c(list(centroids = centroids_sf), extra_cgs),
        annotGeometries = if (length(annots)) annots else NULL,
        spatialGraphs = graphs_arg, unit = u)

    mase_img <- imgData(x)
    if (!is.null(mase_img))
        SpatialExperiment::imgData(sfe) <- mase_img
    sfe
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helpers for splitByCol and aggregate
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @importFrom MultiAssayExperiment experiments
.sfe_to_mase <- function(sfe, template) {
    if (length(experiments(template)) != 1L)
        return(sfe)
    as(sfe, "MultiAssaySpatialExperiment")
}

#' @importFrom MultiAssayExperiment experiments
.sfe_list_to_mase_list <- function(sfe_list, template) {
    if (length(experiments(template)) != 1L)
        return(sfe_list)
    lapply(sfe_list, as, "MultiAssaySpatialExperiment")
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Methods (registered in .onLoad when SFE available)
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.register_SFE_methods <- function() {
    .import_SFE_generics()
    setMethod("dimGeometry", "MultiAssaySpatialExperiment",
        function(x, type = 1L, MARGIN, sample_id = 1L, withDimnames = TRUE) {
            exp <- tryCatch(.materialize_SFE_view(x), error = function(e) e)
            if (inherits(exp, "error"))
                stop("Cannot materialize a SpatialFeatureExperiment view ",
                     "from 'x': ", conditionMessage(exp))
            SpatialFeatureExperiment::dimGeometry(exp,
                type = type, MARGIN = MARGIN, sample_id = sample_id,
                withDimnames = withDimnames)
        },
        where = .sfe_ext_env
    )

    setMethod("dimGeometries", "MultiAssaySpatialExperiment",
        function(x, MARGIN = 2, withDimnames = TRUE) {
            exp <- tryCatch(.materialize_SFE_view(x), error = function(e) e)
            if (inherits(exp, "error"))
                stop("Cannot materialize a SpatialFeatureExperiment view ",
                     "from 'x': ", conditionMessage(exp))
            SpatialFeatureExperiment::dimGeometries(exp,
                MARGIN = MARGIN, withDimnames = withDimnames)
        },
        where = .sfe_ext_env
    )

    setMethod("dimGeometryNames", "MultiAssaySpatialExperiment",
        function(x, MARGIN) {
            exp <- tryCatch(.materialize_SFE_view(x), error = function(e) e)
            if (inherits(exp, "error"))
                stop("Cannot materialize a SpatialFeatureExperiment view ",
                     "from 'x': ", conditionMessage(exp))
            SpatialFeatureExperiment::dimGeometryNames(exp, MARGIN)
        },
        where = .sfe_ext_env
    )

    setMethod("annotGeometry", "MultiAssaySpatialExperiment",
        function(x, type = 1L, sample_id = 1L) {
            exp <- tryCatch(.materialize_SFE_view(x), error = function(e) e)
            if (inherits(exp, "error"))
                stop("Cannot materialize a SpatialFeatureExperiment view ",
                     "from 'x': ", conditionMessage(exp))
            SpatialFeatureExperiment::annotGeometry(exp,
                type = type, sample_id = sample_id)
        },
        where = .sfe_ext_env
    )

    setMethod("annotGeometries", "MultiAssaySpatialExperiment",
        function(x) {
            exp <- tryCatch(.materialize_SFE_view(x), error = function(e) e)
            if (inherits(exp, "error"))
                stop("Cannot materialize a SpatialFeatureExperiment view ",
                     "from 'x': ", conditionMessage(exp))
            SpatialFeatureExperiment::annotGeometries(exp)
        },
        where = .sfe_ext_env
    )

    setMethod("annotGeometryNames", "MultiAssaySpatialExperiment",
        function(x) {
            exp <- tryCatch(.materialize_SFE_view(x), error = function(e) e)
            if (inherits(exp, "error"))
                stop("Cannot materialize a SpatialFeatureExperiment view ",
                     "from 'x': ", conditionMessage(exp))
            SpatialFeatureExperiment::annotGeometryNames(exp)
        },
        where = .sfe_ext_env
    )

    # spatialGraphs/spatialGraph/spatialGraphNames are MASE-native (read
    # colPair()s on the stored experiment directly) rather than routed
    # through .materialize_SFE_view(), and work whether or not SFE is even
    # installed. Only col-margin (MARGIN = 2) graphs are supported -- see
    # the plan doc for why row/annotation margins are out of scope.
    setMethod("spatialGraphs", "MultiAssaySpatialExperiment",
        function(x, MARGIN = NULL, sample_id = "all", name = "all") {
            exps <- experiments(x)
            if (length(exps) != 1L)
                stop("'spatialGraphs' requires exactly one assay; found ",
                     length(exps))
            if (!is.null(MARGIN) && !isTRUE(MARGIN == 2))
                stop("Only col-margin (MARGIN = 2) spatial graphs are ",
                     "currently supported")
            exp <- exps[[1L]]
            cps <- SingleCellExperiment::colPairs(exp)
            nms <- names(cps)
            if (!identical(name, "all"))
                nms <- intersect(nms, name)
            graphs <- stats::setNames(
                lapply(nms, function(nm) .selfhits_to_listw(cps[[nm]], colnames(exp))),
                nms)
            if (length(graphs) == 1L && is.character(name) && !identical(name, "all"))
                graphs[[1L]]
            else
                graphs
        },
        where = .sfe_ext_env
    )

    setMethod("spatialGraph", "MultiAssaySpatialExperiment",
        function(x, type = 1L, MARGIN, sample_id = 1L) {
            exps <- experiments(x)
            if (length(exps) != 1L)
                stop("'spatialGraph' requires exactly one assay; found ",
                     length(exps))
            if (!missing(MARGIN) && !isTRUE(MARGIN == 2))
                stop("Only col-margin (MARGIN = 2) spatial graphs are ",
                     "currently supported")
            exp <- exps[[1L]]
            cps <- SingleCellExperiment::colPairs(exp)
            nm <- if (is.numeric(type)) names(cps)[[type]] else type
            if (is.null(nm) || !nm %in% names(cps))
                stop("No spatial graph named '", type, "' found")
            .selfhits_to_listw(cps[[nm]], colnames(exp))
        },
        where = .sfe_ext_env
    )

    setMethod("spatialGraphNames", "MultiAssaySpatialExperiment",
        function(x, MARGIN, sample_id = 1L) {
            exps <- experiments(x)
            if (length(exps) != 1L)
                return(character(0L))
            if (!missing(MARGIN) && !isTRUE(MARGIN == 2))
                return(character(0L))
            names(SingleCellExperiment::colPairs(exps[[1L]]))
        },
        where = .sfe_ext_env
    )

    setMethod("bbox", "MultiAssaySpatialExperiment",
        function(sfe, sample_id = "all", ...) {
            exp <- tryCatch(.materialize_SFE_view(sfe), error = function(e) e)
            if (inherits(exp, "error"))
                stop("Cannot materialize a SpatialFeatureExperiment view ",
                     "from 'sfe': ", conditionMessage(exp))
            SpatialFeatureExperiment::bbox(exp, sample_id = sample_id, ...)
        },
        where = .sfe_ext_env
    )

    setMethod("unit", "MultiAssaySpatialExperiment",
        function(x) {
            exp <- tryCatch(.materialize_SFE_view(x), error = function(e) e)
            if (inherits(exp, "error"))
                return(NULL)
            SpatialFeatureExperiment::unit(exp)
        },
        where = .sfe_ext_env
    )

    setMethod("localResults", "MultiAssaySpatialExperiment",
        function(x, sample_id = "all", name = "all", features = NULL,
                 colGeometryName = NULL, annotGeometryName = NULL,
                 withDimnames = TRUE, swap_rownames = NULL, ...) {
            stop("'localResults' is not supported on a ",
                 "MultiAssaySpatialExperiment: it is not carried over by ",
                 "the SpatialFeatureExperiment coercion (see ?",
                 "MultiAssaySpatialExperiment-SpatialFeatureExperiment).")
        },
        where = .sfe_ext_env
    )

    setMethod("localResultNames", "MultiAssaySpatialExperiment",
        function(x) character(0L),
        where = .sfe_ext_env
    )

    setMethod("splitByCol", c("MultiAssaySpatialExperiment", "sf"),
        function(x, f, sample_id = "all", colGeometryName = 1L, cover = FALSE) {
            exp <- tryCatch(.materialize_SFE_view(x), error = function(e) e)
            if (inherits(exp, "error"))
                stop("Cannot materialize a SpatialFeatureExperiment view ",
                     "from 'x': ", conditionMessage(exp))
            spl <- SpatialFeatureExperiment::splitByCol(exp,
                f = f, sample_id = sample_id,
                colGeometryName = colGeometryName, cover = cover)
            .sfe_list_to_mase_list(spl, x)
        },
        where = .sfe_ext_env
    )

    setMethod("splitByCol", c("MultiAssaySpatialExperiment", "sfc"),
        function(x, f, sample_id = 1L, colGeometryName = 1L, cover = FALSE) {
            exp <- tryCatch(.materialize_SFE_view(x), error = function(e) e)
            if (inherits(exp, "error"))
                stop("Cannot materialize a SpatialFeatureExperiment view ",
                     "from 'x': ", conditionMessage(exp))
            spl <- SpatialFeatureExperiment::splitByCol(exp,
                f = f, sample_id = sample_id,
                colGeometryName = colGeometryName, cover = cover)
            .sfe_list_to_mase_list(spl, x)
        },
        where = .sfe_ext_env
    )

    setMethod("splitByCol", c("MultiAssaySpatialExperiment", "list"),
        function(x, f, sample_id = "all", colGeometryName = 1L, cover = FALSE) {
            exp <- tryCatch(.materialize_SFE_view(x), error = function(e) e)
            if (inherits(exp, "error"))
                stop("Cannot materialize a SpatialFeatureExperiment view ",
                     "from 'x': ", conditionMessage(exp))
            spl <- SpatialFeatureExperiment::splitByCol(exp,
                f = f, sample_id = sample_id,
                colGeometryName = colGeometryName, cover = cover)
            .sfe_list_to_mase_list(spl, x)
        },
        where = .sfe_ext_env
    )

    setMethod("aggregate", "MultiAssaySpatialExperiment",
        function(x, by = NULL, sample_id = "all", ...) {
            exp <- tryCatch(.materialize_SFE_view(x), error = function(e) e)
            if (inherits(exp, "error"))
                stop("Cannot materialize a SpatialFeatureExperiment view ",
                     "from 'x': ", conditionMessage(exp))
            agg <- SpatialFeatureExperiment::aggregate(exp,
                by = by, sample_id = sample_id, ...)
            .sfe_to_mase(agg, x)
        },
        where = .sfe_ext_env
    )

    setMethod("findDebrisCells", "MultiAssaySpatialExperiment",
        function(x, max_cells = 5, distance_cutoff = 50,
                 BNPARAM = NULL, BPPARAM = BiocParallel::SerialParam()) {
            exp <- tryCatch(.materialize_SFE_view(x), error = function(e) e)
            if (inherits(exp, "error"))
                stop("Cannot materialize a SpatialFeatureExperiment view ",
                     "from 'x': ", conditionMessage(exp))
            SpatialFeatureExperiment::findDebrisCells(exp,
                max_cells = max_cells, distance_cutoff = distance_cutoff,
                BNPARAM = BNPARAM, BPPARAM = BPPARAM)
        },
        where = .sfe_ext_env
    )
}
