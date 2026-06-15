### =========================================================================
### Read 10x Genomics Visium HD data
### -------------------------------------------------------------------------

#' Read 10x Genomics Visium HD Spatial Data
#'
#' @description
#' Load a 10x Genomics Visium HD Space Ranger output directory into a
#' \linkS4class{MultiAssaySpatialExperiment}.
#'
#' @details
#' Visium HD provides binned expression at 2 µm, 8 µm, and 16 µm resolution
#' (codes \code{"002"}, \code{"008"}, \code{"016"}). Each requested bin size
#' becomes a separate assay in the returned object.
#'
#' @param data_dir Character. Path to Space Ranger output directory.
#' @param sample_id Character. Optional sample identifier (default: directory name).
#' @param bin_size Character vector. Bin sizes to load: "002", "008", "016" (default: "008").
#' @param load_boundaries Logical. Whether to load cell segmentation boundaries (default: FALSE).
#' @param images Logical. Whether to load image metadata (default: TRUE).
#' @param load_images Logical. Whether to load image data (default: FALSE).
#' @param min_area Numeric. Minimum polygon area for filtering (default: NULL).
#'
#' @return A \linkS4class{MultiAssaySpatialExperiment} with one experiment per
#'   requested bin size.
#'
#' @seealso \code{\link{readVisiumMASE}}, \code{\link{readXeniumMASE}}
#'
#' @examples
#' \dontrun{
#' ## Requires a Visium HD Space Ranger output directory
#' mase <- readVisiumHDMASE("path/to/visium_hd", bin_size = "008")
#' names(mase)
#' }
#'
#' @importFrom S4Vectors isTRUEorFALSE wmsg
#'
#' @export
readVisiumHDMASE <-
function(data_dir, sample_id = NULL, bin_size = c("008", "002", "016"),
         load_boundaries = FALSE, images = TRUE, load_images = FALSE,
         min_area = NULL)
{    
    ## Argument validation
    data_dir <- normarg_data_dir(data_dir)
    sample_id <- normarg_sample_id(sample_id, data_dir)
    bin_size <- match.arg(bin_size, several.ok = TRUE)
    img_args <- normarg_images(images, load_images)
    min_area <- normarg_min_area(min_area)

    if (!isTRUEorFALSE(load_boundaries))
        stop(wmsg("'load_boundaries' must be TRUE or FALSE"))

    ## File discovery (technology-specific)
    files <- .check_visiumhd_files(data_dir, bin_size)

    ## === COMPONENT READERS ===

    ## For each bin size, read counts and positions
    exp_data <- lapply(bin_size, function(bs) {
        bin_files <- files$bins[[bs]]

        ## 1. Read counts using component reader
        counts <- readHDF5ForMASE(bin_files$matrix, type = "10x")

        ## 2. Read positions using component reader (Parquet!)
        positions <- .read_coldata_parquet(bin_files$positions,
                                           technology = "VisiumHD")

        list(counts = counts, positions = positions)
    })
    names(exp_data) <- paste0("bin_", bin_size)

    ## 3. Read cell boundaries if requested (GeoJSON)
    boundaries <- if (load_boundaries && !is.null(files$boundaries)) {
        .read_geometries_geojson(files$boundaries,  min_area = min_area)
    } else {
        NULL
    }

    ## 4. Read scale factors using component reader
    scale <- if (!is.null(files$scalefactors)) {
        .read_visium_scalefactors(files$scalefactors)
    } else {
        NULL
    }

    ## 5. Read images using component reader
    imgData <- if (img_args$images && !is.null(files$spatial_dir)) {
        .read_visium_images(files$spatial_dir, sample_id,
                            img_args$load_images, scale)
    } else {
        NULL
    }

    ## 6. Assemble multi-bin MASE (technology-specific)
    .build_mase_from_visiumhd(exp_data, boundaries, imgData, sample_id)
}

#' @importFrom S4Vectors wmsg
.check_visiumhd_files <- function(data_dir, bin_sizes) {    
    ## Check for binned matrices
    bins <- list()
    for (bs in bin_sizes) {
        bin_dir <- file.path(data_dir, paste0("binned_outputs/square_", bs, "um"))

        if (!dir.exists(bin_dir)) {
            stop(wmsg("Bin directory not found: ", bin_dir))
        }

        ## Matrix file (HDF5 or 10x directory)
        matrix_file <- file.path(bin_dir, "filtered_feature_bc_matrix.h5")
        if (!file.exists(matrix_file)) {
            matrix_dir <- file.path(bin_dir, "filtered_feature_bc_matrix")
            if (dir.exists(matrix_dir)) {
                matrix_file <- matrix_dir
            } else {
                stop(wmsg("Matrix file not found: ", matrix_file))
            }
        }

        ## Positions file (Parquet)
        positions_file <- file.path(bin_dir, "spatial/tissue_positions.parquet")
        if (!file.exists(positions_file)) {
            stop(wmsg("Positions file not found: ", positions_file))
        }

        bins[[bs]] <- list(matrix = matrix_file, positions = positions_file)
    }

    ## Optional: Cell boundaries (GeoJSON)
    boundaries_file <- file.path(data_dir, "spatial/tissue_positions.geojson")
    if (!file.exists(boundaries_file)) {
        boundaries_file <- NULL
    }

    ## Optional: Scale factors
    scalefactors_file <- file.path(data_dir, "spatial/scalefactors_json.json")
    if (!file.exists(scalefactors_file)) {
        scalefactors_file <- NULL
    }

    ## Optional: Spatial directory for images
    spatial_dir <- file.path(data_dir, "spatial")
    if (!dir.exists(spatial_dir)) {
        spatial_dir <- NULL
    }

    list(bins = bins,
         boundaries = boundaries_file,
         scalefactors = scalefactors_file,
         spatial_dir = spatial_dir)
}

#' @importFrom S4Vectors DataFrame
#' @importFrom MultiAssayExperiment ExperimentList
.build_mase_from_visiumhd <-
function(exp_data, boundaries, imgData, sample_id, labels = RasterLayerList())
{
    ## Create ExperimentList from multiple bin sizes
    exp_list <- lapply(names(exp_data), function(bin_name) {
        exp_data[[bin_name]]$counts
    })
    names(exp_list) <- names(exp_data)

    assay_names <- names(exp_list)
    for (bin_name in assay_names) {
        counts <- exp_list[[bin_name]]
        barcodes <- .observation_ids(counts)
        if (is(counts, "SummarizedExperiment") &&
                (is.null(colnames(counts)) || !length(colnames(counts))))
            colnames(counts) <- barcodes
        exp_list[[bin_name]] <- counts
    }
    experiments <- do.call(ExperimentList, exp_list)

    all_barcodes <- unique(unlist(lapply(exp_list, colnames)))

    colData <- DataFrame(sample_id = rep(sample_id, length(all_barcodes)),
                         row.names = all_barcodes)

    sampleMap_list <- lapply(assay_names, function(bin_name) {
        barcodes <- colnames(exp_list[[bin_name]])
        DataFrame(assay = factor(bin_name, assay_names),
                  primary = barcodes,
                  colname = barcodes)
    })
    sampleMap <- do.call(rbind, sampleMap_list)

    spatialMap_rows <- list()
    shapes_list <- list()

    for (bin_name in names(exp_data)) {
        positions <- exp_data[[bin_name]]$positions

        if (is.null(boundaries)) {
            bin_um <- as.numeric(sub("bin_", "", bin_name))
            radius <- bin_um / 2

            if (!"instance_id" %in% colnames(positions) &&
                    "cell_id" %in% colnames(positions))
                positions$instance_id <- positions$cell_id

            spot_geoms <- .create_spot_geometries(
                positions,
                x_col = "x_centroid",
                y_col = "y_centroid",
                radius = radius
            )

            shapes_list[[bin_name]] <- spot_geoms

            geom_ids <- if ("instance_id" %in% colnames(spot_geoms)) {
                as.character(spot_geoms$instance_id)
            } else {
                colnames(exp_list[[bin_name]])
            }
            bin_map <- .filterSpatialMapByInstances(
                buildSpatialMap(sampleMap,
                                region = bin_name,
                                element_type = "shapes",
                                assays = bin_name),
                geom_ids)
            spatialMap_rows[[length(spatialMap_rows) + 1L]] <- bin_map
        }
    }

    if (!is.null(boundaries)) {
        shapes_list$cells <- boundaries

        boundary_ids <- as.character(boundaries$instance_id)
        for (bin_name in names(exp_data)) {
            cell_map <- .filterSpatialMapByInstances(
                buildSpatialMap(sampleMap,
                                region = "cells",
                                element_type = "shapes",
                                assays = bin_name),
                boundary_ids)
            spatialMap_rows[[length(spatialMap_rows) + 1L]] <- cell_map
        }
    }

    spatialMap <- if (length(spatialMap_rows) > 0L) {
        do.call(.rbindSpatialMaps, spatialMap_rows)
    } else {
        NULL
    }

    shapes <- if (length(shapes_list) > 0L) {
        do.call(ShapesLayerList, shapes_list)
    } else {
        ShapesLayerList()
    }

    MultiAssaySpatialExperiment(experiments = experiments,
                                colData = colData,
                                sampleMap = sampleMap,
                                shapes = shapes,
                                images = .spatialImages_from_imgData(imgData),
                                labels = labels,
                                imgData = imgData,
                                spatialMap = spatialMap)
}
