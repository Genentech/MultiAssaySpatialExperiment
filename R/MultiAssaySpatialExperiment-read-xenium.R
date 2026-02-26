### =========================================================================
### Read 10x Genomics Xenium data
### -------------------------------------------------------------------------

#' Read 10x Genomics Xenium Spatial Data
#'
#' @param data_dir Character. Path to Xenium output directory.
#' @param sample_id Character. Optional sample identifier (default: directory name).
#' @param segmentations Character. Segmentation type: "cell", "nucleus", or "both" (default: "cell").
#' @param images Logical. Whether to load image metadata (default: TRUE).
#' @param load_images Logical. Whether to load image data (default: FALSE).
#' @param add_transcripts Logical. Whether to include transcript data (default: FALSE).
#' @param min_area Numeric. Minimum polygon area for filtering (default: NULL).
#' @param min_phred Integer. Minimum phred score for transcripts (default: 20).
#' @param block_size Numeric. Block size for chunked reading (default: 100MB).
#'
#' @return A MultiAssaySpatialExperiment object.
#'
#' @importFrom S4Vectors isTRUEorFALSE wmsg
#' @importFrom tools file_ext
#'
#' @export
readXeniumMASE <-
function(data_dir, sample_id = NULL,
         segmentations = c("cell", "nucleus", "both"), images = TRUE,
         load_images = FALSE, add_transcripts = FALSE, min_area = NULL,
         min_phred = 20, block_size = NULL)
{
    ## Argument validation
    data_dir <- normarg_data_dir(data_dir)
    sample_id <- normarg_sample_id(sample_id, data_dir)
    segmentations <- match.arg(segmentations)
    img_args <- normarg_images(images, load_images)
    min_area <- normarg_min_area(min_area)
    min_phred <- normarg_min_phred(min_phred)
    block_size <- normarg_block_size(block_size)

    if (!isTRUEorFALSE(add_transcripts))
        stop(wmsg("'add_transcripts' must be TRUE or FALSE"))

    ## File discovery (technology-specific)
    version <- .get_xenium_version(data_dir)
    message(sprintf("Detected Xenium version: %s", version))

    config <- get_technology_config("Xenium")
    files <- .check_xenium_files(data_dir, config)

    ## === COMPONENT READERS ===

    ## 1. Read counts using component reader
    counts <- readHDF5ForMASE(files$matrix)

    ## 2. Read cell metadata using component reader
    ext <- tolower(file_ext(files$cells))
    cells <- if (ext == "parquet") {
        .read_coldata_parquet(files$cells, technology = "Xenium")
    } else {
        .read_coldata_csv(files$cells, technology = "Xenium")
    }

    ## 3. Read cell boundaries using component reader with auto-dispatch
    cell_polys <- if (segmentations %in% c("cell", "both") && !is.null(files$cell_boundaries)) {
        .read_geometries_auto(files$cell_boundaries, layer_name = "cells",
                              min_area = min_area)
    } else NULL
    
    ## 4. Read nucleus boundaries using component reader with auto-dispatch
    nuc_polys <- if (segmentations %in% c("nucleus", "both") && !is.null(files$nucleus_boundaries)) {
        .read_geometries_auto(files$nucleus_boundaries, layer_name = "nuclei",
                              min_area = min_area)
    } else NULL

    ## 5. Read transcripts if requested (component reader)
    transcripts <- if (add_transcripts && !is.null(files$transcripts)) {
        ext_tx <- tolower(file_ext(files$transcripts))
        if (ext_tx == "parquet") {
            tx <- .read_transcripts_parquet(files$transcripts,
                                            x_col = c("x_location", "x"),
                                            y_col = c("y_location", "y"),
                                            feature_col = c("feature_name", "gene"),
                                            qc_col = c("qv", "quality"))
            ## Filter by phred score if requested
            if (!is.null(min_phred) && "qc_score" %in% colnames(tx)) {
                tx <- tx[tx$qc_score >= min_phred, , drop = FALSE]
            }
            tx
        } else if (ext_tx == "csv") {
            .read_transcripts_csv(files$transcripts,
                                  x_col = c("x", "x_location"),
                                  y_col = c("y", "y_location"),
                                  feature_col = c("gene", "feature_name"))
        } else {
            NULL
        }
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
    .build_mase_from_xenium(counts, cells, cell_polys, nuc_polys, transcripts,
                            imgData, sample_id, segmentations)
}

.get_xenium_version <- function(data_dir) {
    exp_file <- file.path(data_dir, "experiment.xenium")
    if (file.exists(exp_file)) {
        return("v2")
    }
    "v1"
}

#' @importFrom S4Vectors wmsg
.check_xenium_files <- function(data_dir, config) {
    matrix_file <- NULL
    for (pattern in config@file_patterns$matrix) {
        candidate <- file.path(data_dir, pattern)
        if (file.exists(candidate) || dir.exists(candidate)) {
            matrix_file <- candidate
            break
        }
    }

    if (is.null(matrix_file)) {
        stop(wmsg(
            "No matrix file found in: ", data_dir, ". ",
            "Tried: ", paste(config@file_patterns$matrix, collapse = ", ")
        ))
    }

    cells_file <- NULL
    for (pattern in config@file_patterns$cells) {
        candidate <- file.path(data_dir, pattern)
        if (file.exists(candidate)) {
            cells_file <- candidate
            break
        }
    }

    if (is.null(cells_file)) {
        stop(wmsg("No cells file found in: ", data_dir, ". ",
                  "Tried: ", paste(config@file_patterns$cells, collapse = ", ")))
    }

    cell_boundaries_file <- NULL
    for (pattern in config@file_patterns$cell_boundaries) {
        candidate <- file.path(data_dir, pattern)
        if (file.exists(candidate)) {
            cell_boundaries_file <- candidate
            break
        }
    }

    nucleus_boundaries_file <- NULL
    for (pattern in config@file_patterns$nucleus_boundaries) {
        candidate <- file.path(data_dir, pattern)
        if (file.exists(candidate)) {
            nucleus_boundaries_file <- candidate
            break
        }
    }

    ## Transcripts file (optional)
    transcripts_file <- NULL
    tx_patterns <- c("transcripts.parquet", "transcripts.csv.gz", "transcripts.csv")
    for (pattern in tx_patterns) {
        candidate <- file.path(data_dir, pattern)
        if (file.exists(candidate)) {
            transcripts_file <- candidate
            break
        }
    }

    image_files <- list.files(data_dir, pattern = config@file_patterns$images,
                              full.names = TRUE)

    list(matrix = matrix_file,
         cells = cells_file,
         cell_boundaries = cell_boundaries_file,
         nucleus_boundaries = nucleus_boundaries_file,
         transcripts = transcripts_file,
         images = image_files)
}

## Technology-specific helpers retained below

#' @importFrom S4Vectors DataFrame
#' @importFrom MultiAssayExperiment ExperimentList
.build_mase_from_xenium <-
function(counts, cells, cell_polys, nuc_polys,transcripts, imgData, sample_id,
         segmentations)
{
    cell_ids <- colnames(counts)

    experiments <- ExperimentList(xenium = counts)

    colData <- DataFrame(sample_id = sample_id, row.names = cell_ids)
    
    sampleMap <- DataFrame(assay = "xenium",
                           primary = cell_ids,
                           colname = cell_ids)

    spatialMap_rows <- list()
    shapes_list <- list()
    points_list <- list()

    if (segmentations %in% c("cell", "both") && !is.null(cell_polys)) {
        cell_ids_in_geom <- as.character(cell_polys$instance_id)
        cell_map <- DataFrame(assay = "xenium",
                              colname = cell_ids,
                              element_type = "shapes",
                              region = "cells",
                              instance_id = cell_ids)
        cell_map <- cell_map[cell_map$instance_id %in% cell_ids_in_geom, ,
                             drop = FALSE]
        spatialMap_rows[[length(spatialMap_rows) + 1L]] <- cell_map
        shapes_list$cells <- cell_polys
    }

    if (segmentations %in% c("nucleus", "both") && !is.null(nuc_polys)) {
        nuc_ids_in_geom <- as.character(nuc_polys$instance_id)
        nuc_map <- DataFrame(assay = "xenium",
                             colname = cell_ids,
                             element_type = "shapes",
                             region = "nuclei",
                             instance_id = cell_ids)
        nuc_map <- nuc_map[nuc_map$instance_id %in% nuc_ids_in_geom, ,
                           drop = FALSE]
        spatialMap_rows[[length(spatialMap_rows) + 1L]] <- nuc_map
        shapes_list$nuclei <- nuc_polys
    }

    ## Add transcripts to points if available
    if (!is.null(transcripts)) {
        ## Create spatialMap entry for transcripts
        tx_map <- DataFrame(assay = "xenium",
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
