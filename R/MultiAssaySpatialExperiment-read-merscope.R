### =========================================================================
### Read Vizgen MERSCOPE data
### -------------------------------------------------------------------------

#' Read Vizgen MERSCOPE Spatial Data
#'
#' @param data_dir Character. Path to MERSCOPE output directory.
#' @param sample_id Character. Optional sample identifier (default: directory name).
#' @param fov_ids Character vector or NULL. FOV IDs to load (default: NULL loads all).
#' @param segmentation Character. Segmentation method: "cellpose", "watershed", or "both" (default: "cellpose").
#' @param load_transcripts Logical. Whether to load transcript data (default: FALSE).
#' @param images Logical. Whether to load image metadata (default: TRUE).
#' @param load_images Logical. Whether to load image data (default: FALSE).
#' @param min_area Numeric. Minimum polygon area for filtering (default: NULL).
#' @param min_qv Numeric. Minimum quality value for transcripts (default: 20).
#'
#' @return A MultiAssaySpatialExperiment object.
#'
#' @importFrom S4Vectors isTRUEorFALSE wmsg
#' @importFrom tools file_ext
#'
#' @export
readMERSCOPEMASE <-
function(data_dir, sample_id = NULL, fov_ids = NULL,
         segmentation = c("cellpose", "watershed", "both"),
         load_transcripts = FALSE, images = TRUE, load_images = FALSE,
         min_area = NULL, min_qv = 20)
{
    ## Argument validation
    data_dir <- normarg_data_dir(data_dir)
    sample_id <- normarg_sample_id(sample_id, data_dir)
    segmentation <- match.arg(segmentation)
    img_args <- normarg_images(images, load_images)
    min_area <- normarg_min_area(min_area)

    if (!isTRUEorFALSE(load_transcripts))
        stop(wmsg("'load_transcripts' must be TRUE or FALSE"))

    ## File discovery (technology-specific)
    files <- .check_merscope_files(data_dir, fov_ids, segmentation)

    message(sprintf("Found %d FOVs", length(files$fovs)))

    ## === COMPONENT READERS ===

    ## 1. Read cell metadata using component reader
    cell_meta <- if (!is.null(files$cell_metadata)) {
        ext <- tolower(tools::file_ext(files$cell_metadata))
        if (ext == "parquet") {
            .read_coldata_parquet(files$cell_metadata, technology = "MERSCOPE")
        } else {
            .read_coldata_csv(files$cell_metadata, technology = "MERSCOPE")
        }
    } else {
        NULL
    }

    ## 2. Read cell boundaries using component reader with auto-dispatch
    boundaries_list <- list()

    if (segmentation %in% c("cellpose", "both") && !is.null(files$cellpose_boundaries)) {
        boundaries_list$cellpose <- .read_geometries_auto(files$cellpose_boundaries,
                                                          layer_name = "cellpose",
                                                          min_area = min_area)
    }

    if (segmentation %in% c("watershed", "both") && !is.null(files$watershed_boundaries)) {
        boundaries_list$watershed <- .read_geometries_auto(files$watershed_boundaries,
                                                           layer_name = "watershed",
                                                           min_area = min_area)
    }

    ## 3. Read transcripts if requested (component reader)
    transcripts <- if (load_transcripts && !is.null(files$transcripts)) {
        message("Reading transcripts (this may take a while for large files)...")
        ext <- tolower(file_ext(files$transcripts))
        tx <- if (ext == "parquet") {
            .read_transcripts_parquet(files$transcripts,
                                      x_col = c("global_x", "x"),
                                      y_col = c("global_y", "y"),
                                      feature_col = c("gene", "feature_name"),
                                      qc_col = c("qv", "quality"))
        } else {
            .read_transcripts_csv(files$transcripts,
                                  x_col = c("global_x", "x"),
                                  y_col = c("global_y", "y"),
                                  feature_col = c("gene", "feature_name"))
        }

        ## Filter by quality if QC column present
        if (!is.null(min_qv) && "qc_score" %in% colnames(tx)) {
            tx <- tx[tx$qc_score >= min_qv, , drop = FALSE]
        }
        tx
    } else {
        NULL
    }

    ## 4. Create synthetic counts from cell metadata (MERSCOPE doesn't have matrix by default)
    counts <- if (!is.null(cell_meta)) {
        .build_merscope_counts(cell_meta, transcripts)
    } else {
        NULL
    }

    ## 5. Read images using component reader
    imgData <- if (img_args$images && length(files$images) > 0L) {
        .read_images_metadata(files$images, sample_id = sample_id,
                              load_images = img_args$load_images)
    } else {
        NULL
    }

    ## 6. Assemble MASE (technology-specific)
    .build_mase_from_merscope(counts, cell_meta, boundaries_list, transcripts,
                              imgData, sample_id, segmentation)
}

.check_merscope_files <- function(data_dir, fov_ids, segmentation) {
    ## Cell metadata file (CSV or Parquet)
    cell_meta_patterns <- c("cell_metadata.csv", "cell_metadata.parquet",
                            "cells.csv", "cells.parquet")
    cell_metadata <- NULL
    for (pattern in cell_meta_patterns) {
        candidate <- file.path(data_dir, "region_0", pattern)
        if (file.exists(candidate)) {
            cell_metadata <- candidate
            break
        }
    }

    ## Cell boundaries files
    cellpose_boundaries <- NULL
    cellpose_patterns <- c("cellpose_micron_space.parquet", 
                           "cell_boundaries_cellpose.parquet",
                           "cell_boundaries.parquet")
    for (pattern in cellpose_patterns) {
        candidate <- file.path(data_dir, "region_0", pattern)
        if (file.exists(candidate)) {
            cellpose_boundaries <- candidate
            break
        }
    }

    watershed_boundaries <- NULL
    watershed_patterns <- c("watershed_micron_space.parquet",
                            "cell_boundaries_watershed.parquet")
    for (pattern in watershed_patterns) {
        candidate <- file.path(data_dir, "region_0", pattern)
        if (file.exists(candidate)) {
            watershed_boundaries <- candidate
            break
        }
    }

    ## Transcript files
    transcripts <- NULL
    tx_patterns <- c("detected_transcripts.csv", "transcripts.parquet",
                     "detected_transcripts.parquet")
    for (pattern in tx_patterns) {
        candidate <- file.path(data_dir, "region_0", pattern)
        if (file.exists(candidate)) {
            transcripts <- candidate
            break
        }
    }

    ## FOV information
    fov_dir <- file.path(data_dir, "region_0", "images")
    fovs <- if (dir.exists(fov_dir)) {
        list.files(fov_dir, pattern = "^mosaic_", full.names = TRUE)
    } else {
        character(0L)
    }

    ## Image files
    image_files <- if (dir.exists(fov_dir)) {
        list.files(fov_dir, pattern = "\\.(tif|tiff|png)$", 
                  full.names = TRUE, recursive = TRUE)
    } else {
        character(0L)
    }

    list(cell_metadata = cell_metadata,
         cellpose_boundaries = cellpose_boundaries,
         watershed_boundaries = watershed_boundaries,
         transcripts = transcripts,
         fovs = fovs,
         images = image_files)
}

#' @importFrom SummarizedExperiment SummarizedExperiment
.build_merscope_counts <- function(cell_meta, transcripts) {
    ## If transcripts available, build count matrix
    if (!is.null(transcripts) && "feature" %in% colnames(transcripts) &&
        "cell_id" %in% colnames(transcripts)) {

        ## Aggregate transcripts by cell and feature
        tx_counts <- table(transcripts$cell_id, transcripts$feature)

        ## Convert to SummarizedExperiment
        se <- SummarizedExperiment(assays = list(counts = tx_counts))
    } else if (!is.null(cell_meta)) {
        ## Create placeholder with cell IDs
        cell_ids <- if ("cell_id" %in% colnames(cell_meta)) {
            cell_meta$cell_id
        } else {
            rownames(cell_meta)
        }

        ## Empty matrix (no counts without transcripts)
        mat <- matrix(0L, nrow = 0L, ncol = length(cell_ids),
                      dimnames = list(NULL, cell_ids))
        se <- SummarizedExperiment(assays = list(counts = mat))
    } else {
        NULL
    }

    se
}

#' @importFrom MultiAssayExperiment ExperimentList
#' @importFrom S4Vectors DataFrame
.build_mase_from_merscope <-
function(counts, cell_meta, boundaries_list, transcripts, imgData, sample_id,
         segmentation)
{
    cell_ids <- if (!is.null(counts)) {
        colnames(counts)
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
        ExperimentList(merscope = counts)
    } else {
        ExperimentList()
    }

    colData <- DataFrame(sample_id = sample_id, row.names = cell_ids)

    sampleMap <- if (length(experiments) > 0L) {
        DataFrame(assay = "merscope", primary = cell_ids, colname = cell_ids)
    } else {
        DataFrame(assay = character(0L), primary = character(0L), colname = character(0L))
    }

    spatialMap_rows <- list()
    shapes_list <- list()
    points_list <- list()

    ## Add boundaries to shapes
    for (seg_name in names(boundaries_list)) {
        boundaries <- boundaries_list[[seg_name]]
        boundary_ids <- as.character(boundaries$instance_id)

        spatialMap_rows[[length(spatialMap_rows) + 1L]] <-
            DataFrame(assay = "merscope",
                      colname = cell_ids,
                      element_type = "shapes",
                      region = seg_name,
                      instance_id = boundary_ids)

        shapes_list[[seg_name]] <- boundaries
    }

    ## Add transcripts to points if available
    if (!is.null(transcripts)) {
        tx_map <- DataFrame(assay = "merscope",
                            colname = NA_character_,
                            element_type = "points",
                            region = "transcripts",
                            instance_id = seq_len(nrow(transcripts)))
        spatialMap_rows[[length(spatialMap_rows) + 1L]] <- tx_map
        points_list$transcripts <- transcripts
    }

    spatialMap <- if (length(spatialMap_rows) > 0L) {
        do.call(rbind, spatialMap_rows)
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

    MultiAssaySpatialExperiment(experiments = experiments,
                                colData = colData,
                                sampleMap = sampleMap,
                                shapes = shapes,
                                points = points,
                                imgData = imgData,
                                spatialMap = spatialMap)
}
