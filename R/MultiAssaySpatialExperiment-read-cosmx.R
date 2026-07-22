### =========================================================================
### Read NanoString CosMx SMI data
### -------------------------------------------------------------------------

#' Read NanoString CosMx SMI Spatial Data
#'
#' @description
#' Load a NanoString CosMx output directory into a
#' \linkS4class{MultiAssaySpatialExperiment}.
#'
#' @details
#' Supports multi-FOV datasets. Use \code{fov_ids} to load specific fields of
#' view, or \code{NULL} to load all FOVs into one object. Transcript coordinates
#' are optional (\code{load_transcripts}).
#'
#' @param data_dir Character. Path to CosMx output directory.
#' @param sample_id Character. Optional sample identifier (default: directory name).
#' @param fov_ids Character vector or NULL. FOV IDs to load (default: NULL loads all).
#' @param load_transcripts Logical. Whether to load transcript data (default: FALSE).
#' @param images Logical. Whether to load image metadata (default: TRUE).
#' @param load_images Logical. Whether to load image data (default: FALSE).
#' @param min_area Numeric. Minimum polygon area for filtering (default: NULL).
#'
#' @return A \linkS4class{MultiAssaySpatialExperiment} object.
#'
#' @seealso \code{\link{readMERSCOPEMASE}}, \code{\link{readXeniumMASE}}
#'
#' @examples
#' ## A minimal mock CosMx output directory is bundled for a runnable example;
#' ## in practice, point 'data_dir' at a real CosMx output directory.
#' dir <- system.file("extdata", "cosmx_mock",
#'                    package = "MultiAssaySpatialExperiment")
#' mase <- readCosMxMASE(dir, images = FALSE, load_transcripts = FALSE)
#' names(mase)
#'
#' @importFrom S4Vectors isTRUEorFALSE wmsg
#' @importFrom SummarizedExperiment colData<- SummarizedExperiment
#'
#' @export
readCosMxMASE <-
function(data_dir, sample_id = NULL, fov_ids = NULL, load_transcripts = FALSE,
         images = TRUE, load_images = FALSE, min_area = NULL)
{
    ## Argument validation
    data_dir <- normarg_data_dir(data_dir)
    sample_id <- normarg_sample_id(sample_id, data_dir)
    img_args <- normarg_images(images, load_images)
    min_area <- normarg_min_area(min_area)

    if (!isTRUEorFALSE(load_transcripts))
        stop(wmsg("'load_transcripts' must be TRUE or FALSE"))

    ## File discovery (technology-specific)
    files <- .check_cosmx_files(data_dir, fov_ids)

    message(sprintf("Found %d FOVs", length(files$fovs)))

    ## === COMPONENT READERS ===

    ## 1. Read cell metadata using component reader
    cell_meta <- if (!is.null(files$cell_metadata)) {
        .read_coldata_csv(files$cell_metadata, technology = "CosMx")
    } else {
        NULL
    }

    ## 2. Read cell boundaries using component reader (GeoParquet!)
    boundaries <- if (!is.null(files$boundaries)) {
        ## CosMx typically provides GeoParquet boundaries
        .read_geometries_auto(files$boundaries, layer_name = "cells",
                              min_area = min_area)
    } else {
        NULL
    }

    ## 3. Read expression matrix (counts)
    counts <- if (!is.null(files$expression)) {
        ## Read CSV and build count matrix
        df <- readCSVForMASE(files$expression)

        ## Extract gene names from first column
        gene_names <- df[[1L]]
        df[[1L]] <- NULL

        ## Convert to matrix
        mat <- as.matrix(df)
        rownames(mat) <- gene_names

        ## Create SummarizedExperiment
        se <- SummarizedExperiment(assays = list(counts = mat))
        colData(se)$sample_id <- sample_id

        se
    } else {
        NULL
    }

    ## 4. Read transcripts if requested (component reader)
    transcripts <- if (load_transcripts && !is.null(files$transcripts)) {
        message("Reading transcripts (this may take a while for large files)...")
        ext <- tolower(tools::file_ext(files$transcripts))
        if (ext == "parquet") {
            .read_transcripts_parquet(files$transcripts,
                                      x_col = c("x_global_px", "x"),
                                      y_col = c("y_global_px", "y"),
                                      feature_col = c("target", "gene"),
                                      cell_col = c("cell_ID", "cell_id"))
        } else {
            .read_transcripts_csv(files$transcripts,
                                  x_col = c("x_global_px", "x"),
                                  y_col = c("y_global_px", "y"),
                                  feature_col = c("target", "gene"))
        }
    } else {
        NULL
    }

    ## 5. Read FOV positions metadata
    fov_meta <- if (!is.null(files$fov_positions)) {
        .read_metadata_csv(files$fov_positions)
    } else {
        NULL
    }

    ## 6. Read images using component reader
    imgData <- if (img_args$images && length(files$images) > 0L) {
        .read_images_metadata(files$images,
                              sample_id = sample_id,
                              load_images = img_args$load_images)
    } else {
        NULL
    }

    ## 7. Assemble MASE (technology-specific)
    .build_mase_from_cosmx(counts, cell_meta, boundaries, 
                           transcripts, fov_meta, imgData, sample_id)
}

#' @importFrom utils read.csv
.check_cosmx_files <- function(data_dir, fov_ids) {
    ## Cell metadata file
    cell_meta_patterns <- c("*_metadata_file.csv", "*_cell_metadata.csv",
                            "cell_metadata.csv")
    cell_metadata <- NULL
    for (pattern in cell_meta_patterns) {
        candidates <- list.files(data_dir, pattern = pattern, full.names = TRUE,
                                 recursive = TRUE)
        if (length(candidates) > 0L) {
            cell_metadata <- candidates[1L]
            break
        }
    }

    ## Cell boundaries (GeoParquet or Parquet)
    boundaries_patterns <- c("*_cell_boundaries_sf.parquet",
                             "*_boundaries.parquet",
                             "cell_boundaries.parquet",
                             "cell_boundaries.csv",
                             "*_cell_boundaries.csv")
    boundaries <- NULL
    for (pattern in boundaries_patterns) {
        candidates <- list.files(data_dir, pattern = pattern, full.names = TRUE,
                                 recursive = TRUE)
        if (length(candidates) > 0L) {
            boundaries <- candidates[1L]
            break
        }
    }

    ## Expression matrix
    expression_patterns <- c("*_exprMat_file.csv", "*_expression.csv",
                             "expression_matrix.csv")
    expression <- NULL
    for (pattern in expression_patterns) {
        candidates <- list.files(data_dir, pattern = pattern, full.names = TRUE,
                                 recursive = TRUE)
        if (length(candidates) > 0L) {
            expression <- candidates[1L]
            break
        }
    }

    ## Transcripts
    tx_patterns <- c("*_tx_file.csv", "*_transcripts.csv",
                     "transcripts.parquet")
    transcripts <- NULL
    for (pattern in tx_patterns) {
        candidates <- list.files(data_dir, pattern = pattern, full.names = TRUE,
                                 recursive = TRUE)
        if (length(candidates) > 0L) {
            transcripts <- candidates[1L]
            break
        }
    }

    ## FOV positions
    fov_patterns <- c("*_fov_positions_file.csv", "fov_positions.csv")
    fov_positions <- NULL
    for (pattern in fov_patterns) {
        candidates <- list.files(data_dir, pattern = pattern, full.names = TRUE,
                                 recursive = TRUE)
        if (length(candidates) > 0L) {
            fov_positions <- candidates[1L]
            break
        }
    }

    ## FOV list
    fovs <- if (!is.null(fov_positions)) {
        fov_df <- read.csv(fov_positions)
        if ("fov" %in% colnames(fov_df)) {
            unique(fov_df$fov)
        } else {
            character(0L)
        }
    } else {
        character(0L)
    }

    ## Image files
    image_files <- list.files(data_dir, pattern = "\\.(tif|tiff|jpg|jpeg)$",
                              full.names = TRUE, recursive = TRUE)

    list(cell_metadata = cell_metadata,
         boundaries = boundaries,
         expression = expression,
         transcripts = transcripts,
         fov_positions = fov_positions,
         fovs = fovs,
         images = image_files)
}

#' @importFrom MultiAssayExperiment ExperimentList
#' @importFrom S4Vectors DataFrame
.build_mase_from_cosmx <-
function(counts, cell_meta, boundaries, transcripts, fov_meta, imgData, sample_id,
         labels = RasterLayerList()) {
    cell_ids <- if (!is.null(counts)) {
        .observation_ids(counts)
    } else if (!is.null(cell_meta)) {
        if ("cell_id" %in% colnames(cell_meta)) {
            as.character(cell_meta$cell_id)
        } else {
            rownames(cell_meta)
        }
    } else {
        character(0L)
    }

    experiments <- if (!is.null(counts)) {
        ExperimentList(cosmx = counts)
    } else {
        ExperimentList()
    }

    if (is(counts, "SummarizedExperiment") &&
            (is.null(colnames(counts)) || !length(colnames(counts))))
        colnames(counts) <- cell_ids

    colData <- DataFrame(sample_id = rep(sample_id, length(cell_ids)),
                         row.names = cell_ids)

    sampleMap <- if (length(experiments) > 0L) {
        DataFrame(assay = factor("cosmx", "cosmx"),
                  primary = cell_ids,
                  colname = cell_ids)
    } else {
        DataFrame(assay = factor(levels = "cosmx"),
                  primary = character(0L),
                  colname = character(0L))
    }

    spatialMap_rows <- list()
    shapes_list <- list()
    points_list <- list()

    ## Add cell boundaries to shapes
    if (!is.null(boundaries)) {
        boundary_ids <- as.character(boundaries$instance_id)

        cell_map <- .filterSpatialMapByInstances(
            buildSpatialMap(sampleMap, region = "cells", element_type = "shapes"),
            boundary_ids)
        spatialMap_rows[[length(spatialMap_rows) + 1L]] <- cell_map

        shapes_list$cells <- boundaries
    }

    ## Add transcripts to points if available
    if (!is.null(transcripts)) {
        spatialMap_rows[[length(spatialMap_rows) + 1L]] <-
            .transcriptSpatialMapRow("cosmx", nrow(transcripts))
        points_list$transcripts <- transcripts
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

    points <- if (length(points_list) > 0L) {
        do.call(PointsLayerList, points_list)
    } else {
        PointsLayerList()
    }

    ## Add FOV metadata to metadata slot if available
    metadata_list <- list()
    if (!is.null(fov_meta)) {
        metadata_list$fov_positions <- fov_meta
    }

    MultiAssaySpatialExperiment(experiments = experiments,
                                colData = colData,
                                sampleMap = sampleMap,
                                metadata = metadata_list,
                                shapes = shapes,
                                points = points,
                                images = .spatialImages_from_imgData(imgData),
                                labels = labels,
                                imgData = imgData,
                                spatialMap = spatialMap)
}
