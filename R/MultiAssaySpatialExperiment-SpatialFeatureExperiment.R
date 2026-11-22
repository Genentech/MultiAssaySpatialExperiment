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

#' Interoperating with SpatialFeatureExperiment
#'
#' @description
#' When the \pkg{SpatialFeatureExperiment} package is loaded, a
#' \code{MultiAssaySpatialExperiment} converts to and from a
#' \code{SpatialFeatureExperiment} and answers its generics. Nothing is
#' registered until \pkg{SpatialFeatureExperiment} is loaded, and it does not
#' matter which package is loaded first.
#'
#' @details
#' \code{as(sfe, "MultiAssaySpatialExperiment")} stores a plain
#' \code{SingleCellExperiment}, with the spatial content in the layers. The
#' centroids become the \code{coordinates} points layer. Other
#' \code{colGeometries} become shapes layers linked through \code{spatialMap},
#' and \code{annotGeometries} become shapes layers with no \code{spatialMap}
#' rows, which is how they are told apart on the way back. Column graphs become
#' \code{colPair}s; with several samples the graphs are merged into one
#' \code{colPair} and split per sample again by
#' \code{as(mase, "SpatialFeatureExperiment")}.
#'
#' Not carried over, with a warning that names them: row geometries, row
#' graphs, annotation graphs, and local results.
#'
#' @section Methods:
#' \code{dimGeometry}, \code{annotGeometry}, \code{bbox},
#' \code{splitByCol}, \code{aggregate} and \code{findDebrisCells} work by
#' building a \code{SpatialFeatureExperiment} from the object on each call
#' (the same conversion as \code{as()}) and calling the SFE method on it, so a
#' call costs a conversion. \code{aggregate}, \code{splitByCol} and
#' \code{findDebrisCells} convert the result back, so they return a
#' \code{MultiAssaySpatialExperiment} when the object has a single assay.
#' These methods need exactly one assay.
#'
#' \code{spatialGraphs}, \code{spatialGraph} and \code{spatialGraphNames} read
#' the \code{colPair}s, one sample at a time, and return what the
#' \pkg{SpatialFeatureExperiment} methods return: a list by sample (of
#' \code{row}, \code{col} and \code{annot} graphs, or of the column graphs
#' when \code{MARGIN = 2}), and a single \code{listw} from \code{spatialGraph}.
#' Only column graphs exist, so \code{row} and \code{annot} are \code{NULL}.
#' Like the others they are registered once \pkg{SpatialFeatureExperiment} is
#' loaded.
#'
#' \code{localResultNames} is always empty and \code{localResults} returns an
#' empty list, because local results are not carried over; asking for a named
#' result is an error.
#'
#' Two results differ from those of the original \code{SpatialFeatureExperiment}
#' because the conversion rebuilds the object. \code{dimGeometry} returns the
#' geometries with an \code{instance_id} column, where SFE has
#' \code{sample_id}. \code{bbox} is computed from the geometries, so it can
#' differ from a bounding box stored on the original SFE.
#'
#' \code{unit} is the SFE unit when the single assay is a
#' \code{SpatialFeatureExperiment}, else the \code{spatialUnit} in the
#' metadata, else \code{"full_res_image_pixel"}. An SFE with no unit set gives
#' \code{NULL}, so the two can differ for such an object.
#'
#' @return
#' The coercions return the converted object. The methods return what the
#' \pkg{SpatialFeatureExperiment} method returns, with the exceptions above.
#'
#' @author Patrick Aboyoun
#' @name MultiAssaySpatialExperiment-SpatialFeatureExperiment
#' @aliases coerce,SpatialFeatureExperiment,MultiAssaySpatialExperiment-method
#' @aliases coerce,MultiAssaySpatialExperiment,SpatialFeatureExperiment-method
#' @aliases dimGeometry,MultiAssaySpatialExperiment-method
#' @aliases dimGeometries,MultiAssaySpatialExperiment-method
#' @aliases dimGeometryNames,MultiAssaySpatialExperiment-method
#' @aliases annotGeometries,MultiAssaySpatialExperiment-method
#' @aliases annotGeometryNames,MultiAssaySpatialExperiment-method
#' @aliases spatialGraphs,MultiAssaySpatialExperiment-method
#' @aliases spatialGraph,MultiAssaySpatialExperiment-method
#' @aliases spatialGraphNames,MultiAssaySpatialExperiment-method
#' @aliases unit,MultiAssaySpatialExperiment-method
#' @aliases aggregate,MultiAssaySpatialExperiment-method
#' @aliases splitByCol,MultiAssaySpatialExperiment,sf-method
#' @aliases splitByCol,MultiAssaySpatialExperiment,sfc-method
#' @aliases splitByCol,MultiAssaySpatialExperiment,list-method
#' @aliases annotGeometry,MultiAssaySpatialExperiment-method
#' @aliases bbox,MultiAssaySpatialExperiment-method
#' @aliases localResults,MultiAssaySpatialExperiment-method
#' @aliases localResultNames,MultiAssaySpatialExperiment-method
#' @aliases findDebrisCells,MultiAssaySpatialExperiment-method
#' @seealso \link{coerce-SpatialData} for the spatialdataR counterpart.
#' @examples
#' ## Needs SpatialFeatureExperiment to be loaded
#' if (requireNamespace("SpatialFeatureExperiment", quietly = TRUE)) {
#'     library(SpatialFeatureExperiment)
#'     sfe <- as(SpatialExperiment::SpatialExperiment(
#'         list(counts = matrix(1:12, 3, 4,
#'             dimnames = list(letters[1:3], LETTERS[1:4]))),
#'         spatialCoords = cbind(x = 1:4, y = 4:1),
#'         sample_id = "s1"), "SpatialFeatureExperiment")
#'     mase <- as(sfe, "MultiAssaySpatialExperiment")
#'     bbox(mase)
#'     back <- as(mase, "SpatialFeatureExperiment")
#' }
NULL


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helpers: sf <-> MASE-native DataFrame, listw <-> SelfHits
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Older SFE objects (e.g. the SFEData datasets) store the unit as
# "full_res_image_pixels", which the SFE constructor no longer accepts.
.normalize_sfe_unit <- function(u) {
    if (is.null(u)) return(NULL)
    sub("^full_res_image_pixels$", "full_res_image_pixel", u)
}

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
    if (!is.null(cn)) {
        if (!is.null(rownames(sfobj)) && all(cn %in% rownames(sfobj)))
            sfobj <- sfobj[cn, , drop = FALSE]
        else if ("instance_id" %in% colnames(d) && nrow(sfobj) == length(cn))
            warning("the instance ids of the layer do not match the column ",
                    "names of the assay, so its rows are taken in the order ",
                    "stored", call. = FALSE)
    }
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
    if (style %in% c("W", "B", "C", "U", "S", "minmax"))
        return(spdep::nb2listw(nb, glist = wl, style = style, zero.policy = TRUE))
    ## SFE's "raw" style keeps the weights as given, which nb2listw() does not
    ## accept, so the listw is assembled directly, as SFE does.
    listw <- list(style = style, neighbours = nb, weights = wl)
    attr(listw, "region.id") <- attr(nb, "region.id")
    attr(listw, "zero.policy") <- TRUE
    class(listw) <- c("listw", "nb")
    listw
}

# One graph over all columns from the graphs of each sample. 'idx' holds, for
# each sample, the positions of its columns in the whole object.
#' @importFrom S4Vectors SelfHits from to mcols mcols<- metadata metadata<-
.merge_selfhits <- function(hits, idx, nnode, samples = names(hits)) {
    fr <- unlist(Map(function(h, i) i[from(h)], hits, idx), use.names = FALSE)
    to_ <- unlist(Map(function(h, i) i[to(h)], hits, idx), use.names = FALSE)
    w <- unlist(lapply(hits, function(h) mcols(h)[["weight"]]), use.names = FALSE)
    o <- order(fr, to_)
    out <- SelfHits(fr[o], to_[o], nnode = nnode)
    mcols(out) <- DataFrame(weight = w[o])
    ## SFE keeps a weight style per sample; the merged graph keeps them all
    styles <- vapply(hits, function(h) metadata(h)[["style"]] %||% "W", "")
    if (!is.null(samples)) names(styles) <- samples
    metadata(out) <- list(style = styles[[1L]], styles = styles)
    out
}

# The graph among the columns at positions 'idx' of the whole object.
.subset_selfhits <- function(hits, idx, sample = NULL) {
    fr <- match(from(hits), idx)
    to_ <- match(to(hits), idx)
    keep <- !is.na(fr) & !is.na(to_)
    out <- SelfHits(fr[keep], to_[keep], nnode = length(idx))
    mcols(out) <- DataFrame(weight = mcols(hits)[["weight"]][keep])
    md <- metadata(hits)
    styles <- md[["styles"]]
    if (!is.null(sample) && !is.null(styles) && sample %in% names(styles))
        md[["style"]] <- styles[[sample]]
    md[["styles"]] <- NULL
    metadata(out) <- md
    out
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
            sids <- if ("sample_id" %in% colnames(SummarizedExperiment::colData(from)))
                as.character(SummarizedExperiment::colData(from)[["sample_id"]])
            else rep("sample01", n)
            samples <- unique(sids)

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
            sm_pts[["__instance_id__"]] <- pt_df[["instance_id"]]
            spatial_map <- buildSpatialMap(sm_pts, region = "coordinates",
                element_type = "points", instance_from = "__instance_id__")

            # Every other colGeometry -> a per-cell 'shapes' entry, with
            # matching spatialMap rows.
            shapes_list <- list()
            other_cg <- setdiff(SpatialFeatureExperiment::colGeometryNames(from),
                                "centroids")
            for (nm in other_cg) {
                cg <- SpatialFeatureExperiment::colGeometry(from, type = nm,
                    sample_id = "all")
                cg_inst <- rownames(cg)
                if (is.null(cg_inst) || length(cg_inst) != nrow(cg))
                    cg_inst <- cn
                shapes_list[[nm]] <- .sf_to_shapes_df(cg, cg_inst)
                sm_cg <- sm
                sm_cg[["__instance_id__"]] <- cg_inst
                spatial_map <- .rbindSpatialMaps(spatial_map,
                    buildSpatialMap(sm_cg, region = nm, element_type = "shapes",
                        instance_from = "__instance_id__"))
            }

            # annotGeometries -> 'shapes' entries with NO spatialMap rows;
            # that absence is what marks them as annotation-level on the way
            # back (they are not one-row-per-cell).
            ag_names <- SpatialFeatureExperiment::annotGeometryNames(from)
            shapes_list[ag_names] <- lapply(ag_names, function(nm)
                .sf_to_shapes_df(SpatialFeatureExperiment::annotGeometry(from,
                    type = nm, sample_id = "all")))

            # Col-margin spatial graphs -> colPair()s on the stored SCE. SFE keeps
            # one graph per sample, so those are merged into one over all columns.
            idx <- lapply(samples, function(s) which(sids == s))
            graph_names <- unique(unlist(lapply(samples, function(s)
                SpatialFeatureExperiment::colGraphNames(from, sample_id = s))))
            for (nm in graph_names) {
                parts <- lapply(samples, function(s) {
                    if (nm %in% SpatialFeatureExperiment::colGraphNames(from, sample_id = s))
                        .listw_to_selfhits(SpatialFeatureExperiment::colGraph(
                            from, type = nm, sample_id = s))
                })
                have <- !vapply(parts, is.null, NA)
                colPair(sce, nm) <- .merge_selfhits(parts[have], idx[have], n, samples[have])
            }

            # Said out loud, because nothing here carries them over.
            lost <- c(
                if (length(rg <- SpatialFeatureExperiment::rowGeometryNames(from)))
                    paste0("rowGeometries (", toString(rg), ")"),
                if (length(rgr <- unique(unlist(lapply(samples, function(s)
                        SpatialFeatureExperiment::rowGraphNames(from, sample_id = s))))))
                    paste0("row graphs (", toString(rgr), ")"),
                if (length(agr <- unique(unlist(lapply(samples, function(s)
                        SpatialFeatureExperiment::annotGraphNames(from, sample_id = s))))))
                    paste0("annotation graphs (", toString(agr), ")"))
            if (length(lost))
                warning("not carried over by this coercion: ",
                        paste(lost, collapse = "; "), call. = FALSE)

            lr <- tryCatch(SpatialFeatureExperiment::localResults(from),
                error = function(e) NULL)
            if (!is.null(lr) && length(lr))
                message("'localResults' on the SpatialFeatureExperiment are ",
                        "not carried over by this coercion; access them on ",
                        "the original object before coercing.")

            experiments <- ExperimentList(stats::setNames(list(sce), assay_name))
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

            ## the object's own metadata goes with it
            if (length(metadata(from)))
                metadata(mase) <- metadata(from)
            u <- tryCatch(SpatialFeatureExperiment::unit(from),
                error = function(e) NULL)
            if (!is.null(u))
                metadata(mase)[["spatialUnit"]] <- .normalize_sfe_unit(u)
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
#' @importFrom S4Vectors metadata metadata<-
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

    u <- .normalize_sfe_unit(metadata(x)[["spatialUnit"]])
    if (is.null(u)) u <- "full_res_image_pixel"

    sids <- as.character(colData(exp)[["sample_id"]])
    if (length(sids) == 0L)
        sids <- rep("sample01", length(cn))
    sample_ids <- unique(sids)
    ## one graph per sample, over that sample's columns only
    graphs_arg <- if (length(cps))
            stats::setNames(lapply(sample_ids, function(sid) {
                idx <- which(sids == sid)
                col <- stats::setNames(
                    lapply(names(cps), function(nm)
                        .selfhits_to_listw(.subset_selfhits(cps[[nm]], idx, sid), cn[idx])),
                    names(cps))
                list(row = NULL, col = col, annot = NULL)
            }), sample_ids)
        else NULL

    sfe <- SpatialFeatureExperiment::SpatialFeatureExperiment(
        assays = assays(exp), rowData = rowData(exp), colData = colData(exp),
        colGeometries = c(list(centroids = centroids_sf), extra_cgs),
        annotGeometries = if (length(annots)) annots else NULL,
        spatialGraphs = graphs_arg, unit = u)

    md <- metadata(x)
    md[["spatialUnit"]] <- NULL
    if (length(md))
        metadata(sfe) <- md
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
### Spatial graphs: one merged colPair per name, read back one sample at a time
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.sfe_graph_experiment <- function(x, what) {
    exps <- experiments(x)
    if (length(exps) != 1L)
        stop("'", what, "' requires exactly one assay; found ", length(exps),
             call. = FALSE)
    exps[[1L]]
}

.sfe_sample_ids <- function(exp) {
    sids <- as.character(SummarizedExperiment::colData(exp)[["sample_id"]])
    if (length(sids) == 0L) rep("sample01", ncol(exp)) else sids
}

# The sample ids selected by 'sample_id' ("all", names, or positions), as
# SpatialFeatureExperiment reads them.
.sfe_graph_samples <- function(exp, sample_id) {
    have <- unique(.sfe_sample_ids(exp))
    if (identical(sample_id, "all")) return(have)
    sid <- if (is.numeric(sample_id)) have[sample_id] else sample_id
    if (anyNA(sid) || !all(sid %in% have))
        stop("sample_id not found: ", toString(sample_id), call. = FALSE)
    sid
}

.sfe_graph_of <- function(exp, name, sample_id) {
    idx <- which(.sfe_sample_ids(exp) == sample_id)
    .selfhits_to_listw(
        .subset_selfhits(SingleCellExperiment::colPairs(exp)[[name]], idx,
                         sample_id),
        colnames(exp)[idx])
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
        function(x, MARGIN, withDimnames = TRUE) {
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
    # loaded. Only col-margin (MARGIN = 2) graphs are supported -- see
    # the plan doc for why row/annotation margins are out of scope.
    setMethod("spatialGraphs", "MultiAssaySpatialExperiment",
        function(x, MARGIN = NULL, sample_id = "all", name = "all") {
            exp <- .sfe_graph_experiment(x, "spatialGraphs")
            if (!is.null(MARGIN) && !isTRUE(MARGIN == 2))
                stop("Only col-margin (MARGIN = 2) spatial graphs are ",
                     "currently supported")
            if (is.null(MARGIN) && !identical(name, "all"))
                stop("Cannot get graphs of the same name across different ",
                     "margins")
            sids <- .sfe_graph_samples(exp, sample_id)
            nms <- names(SingleCellExperiment::colPairs(exp))
            if (!identical(name, "all")) {
                if (!all(name %in% nms))
                    stop("No spatial graph named ",
                         toString(setdiff(name, nms)), " found")
                nms <- name
            }
            lapply(stats::setNames(sids, sids), function(sid) {
                col <- stats::setNames(
                    lapply(nms, function(nm) .sfe_graph_of(exp, nm, sid)), nms)
                if (is.null(MARGIN)) list(row = NULL, col = col, annot = NULL)
                else col
            })
        },
        where = .sfe_ext_env
    )

    setMethod("spatialGraph", "MultiAssaySpatialExperiment",
        function(x, type = 1L, MARGIN, sample_id = 1L) {
            exp <- .sfe_graph_experiment(x, "spatialGraph")
            if (!missing(MARGIN) && !isTRUE(MARGIN == 2))
                stop("Only col-margin (MARGIN = 2) spatial graphs are ",
                     "currently supported")
            nms <- names(SingleCellExperiment::colPairs(exp))
            nm <- if (is.numeric(type)) nms[type] else type
            if (anyNA(nm) || length(nm) != 1L || !nm %in% nms)
                stop("No spatial graph named '", type, "' found")
            sid <- .sfe_graph_samples(exp, sample_id)
            if (length(sid) != 1L)
                stop("Only one sample can be specified at a time.")
            .sfe_graph_of(exp, nm, sid)
        },
        where = .sfe_ext_env
    )

    setMethod("spatialGraphNames", "MultiAssaySpatialExperiment",
        function(x, MARGIN, sample_id = 1L) {
            exp <- .sfe_graph_experiment(x, "spatialGraphNames")
            if (!missing(MARGIN) && !isTRUE(MARGIN == 2))
                return(character(0L))
            .sfe_graph_samples(exp, sample_id)
            names(SingleCellExperiment::colPairs(exp))
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
            exps <- experiments(x)
            if (length(exps) == 1L && is(exps[[1L]], "SpatialFeatureExperiment"))
                return(SpatialFeatureExperiment::unit(exps[[1L]]))
            .normalize_sfe_unit(metadata(x)[["spatialUnit"]]) %||%
                "full_res_image_pixel"
        },
        where = .sfe_ext_env
    )

    setMethod("localResults", "MultiAssaySpatialExperiment",
        function(x, sample_id = "all", name = "all", features = NULL,
                 colGeometryName = NULL, annotGeometryName = NULL,
                 withDimnames = TRUE, swap_rownames = NULL, ...) {
            ## No local results are carried over, so there is nothing to
            ## return (an empty list, as SFE gives when none are stored);
            ## asking for a specific result cannot be satisfied.
            if (!identical(name, "all") || !is.null(features))
                stop("'localResults' are not carried over by the ",
                     "SpatialFeatureExperiment coercion, so '", name[1L],
                     "' is not available on a MultiAssaySpatialExperiment ",
                     "(see ?MultiAssaySpatialExperiment-",
                     "SpatialFeatureExperiment).")
            list()
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
            res <- SpatialFeatureExperiment::findDebrisCells(exp,
                max_cells = max_cells, distance_cutoff = distance_cutoff,
                BNPARAM = BNPARAM, BPPARAM = BPPARAM)
            .sfe_to_mase(res, x)
        },
        where = .sfe_ext_env
    )
}
