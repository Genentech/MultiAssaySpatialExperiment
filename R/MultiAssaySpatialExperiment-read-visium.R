### =========================================================================
### Read 10x Genomics Visium data
### -------------------------------------------------------------------------

#' Read 10x Genomics Visium Spatial Data
#'
#' @description
#' Load a 10x Genomics Visium Space Ranger output directory into a
#' \linkS4class{MultiAssaySpatialExperiment}.
#'
#' @details
#' Reads spot-by-gene counts, tissue positions, spot geometries, and optional
#' H&E image metadata (\code{images}). Coordinates may be returned in pixels or
#' microns (\code{unit}). For Visium HD data, use \code{\link{readVisiumHDMASE}}.
#'
#' @param data_dir Character. Path to Space Ranger output directory.
#' @param sample_id Character. Optional sample identifier (default: directory name).
#' @param type Character. Matrix format: "HDF5" or "sparse" (default: "HDF5").
#' @param data Character. Data type: "filtered" or "raw" (default: "filtered").
#' @param images Logical. Whether to load image metadata (default: TRUE).
#' @param load_images Logical. Whether to load image data (default: FALSE).
#' @param unit Character. Coordinate units: "pixel" or "micron" (default: "pixel").
#' @param min_area Numeric. Minimum polygon area for filtering (default: NULL).
#' @param block_size Numeric. Block size for chunked reading (default: 100MB).
#'
#' @return A \linkS4class{MultiAssaySpatialExperiment} object.
#'
#' @seealso \code{\link{readVisiumHDMASE}}, \code{\link{readXeniumMASE}}
#'
#' @export
readVisiumMASE <-
function(data_dir, sample_id = NULL, type = c("HDF5", "sparse"),
         data = c("filtered", "raw"), images = TRUE, load_images = FALSE,
         unit = c("pixel", "micron"), min_area = NULL, block_size = NULL)
{
    ## Argument validation
    data_dir <- normarg_data_dir(data_dir)
    sample_id <- normarg_sample_id(sample_id, data_dir)
    type <- normarg_matrix_format(type)
    data <- normarg_data_type(data)
    unit <- normarg_unit(unit)
    img_args <- normarg_images(images, load_images)
    min_area <- normarg_min_area(min_area)
    block_size <- normarg_block_size(block_size)

    config <- get_technology_config("Visium")

    ## File discovery (technology-specific)
    files <- .check_visium_files(data_dir, data, config)

    ## === COMPONENT READERS ===

    ## 1. Read counts using component reader
    counts <- readHDF5ForMASE(files$matrix)

    ## 2. Read positions using component reader
    positions <- .read_coldata_csv(files$positions, technology = "Visium")

    ## 3. Read scale factors using component reader
    scale <- .read_visium_scalefactors(files$scalefactors)

    ## 4. Convert units if requested (technology-specific)
    if (unit == "micron") {
        positions <- .convert_positions_to_micron(positions, scale, config)
    }

    ## 5. Create spot geometries using component reader
    spot_radius <- .get_spot_radius(scale, unit, config)

    ## Prepare positions for geometry creation (add instance_id)
    positions$instance_id <- positions$cell_id

    geometries <- .create_spot_geometries(positions,
                                          x_col = "x_centroid",
                                          y_col = "y_centroid",
                                          radius = spot_radius,
                                          min_area = min_area)

    ## Note: Other columns (in_tissue, array_row, array_col) are preserved by sf conversion

    ## 6. Read images using component reader
    imgData <- if (img_args$images) {
        .read_visium_images(files$spatial_dir, sample_id, img_args$load_images,
                            scale)
    } else {
        NULL
    }

    ## 7. Assemble MASE (technology-specific)
    .build_mase_from_visium(counts, positions, geometries, imgData, sample_id)
}

#' @importFrom S4Vectors wmsg
.check_visium_files <- function(data_dir, data, config) {
    matrix_dir <- file.path(data_dir, paste0(data, "_feature_bc_matrix"))
    if (!dir.exists(matrix_dir)) {
        stop(wmsg("Matrix directory not found: ", matrix_dir, ". ",
                  "Expected a Space Ranger output directory with '",
                  paste0(data, "_feature_bc_matrix"), "' subdirectory."))
    }

    h5_file <- file.path(matrix_dir, "matrix.h5")
    mtx_file <- file.path(matrix_dir, "matrix.mtx.gz")

    has_h5 <- file.exists(h5_file)
    has_mtx <- file.exists(mtx_file)

    if (!has_h5 && !has_mtx) {
        stop(wmsg("No matrix file found in: ", matrix_dir, ". ",
                  "Expected either 'matrix.h5' (HDF5) or 'matrix.mtx.gz' (sparse)."))
    }

    matrix_path <- if (has_h5) h5_file else matrix_dir

    spatial_dir <- file.path(data_dir, "spatial")
    if (!dir.exists(spatial_dir)) {
        stop(wmsg("Spatial directory not found: ", spatial_dir, ". ",
                  "Expected a Space Ranger output directory with 'spatial/' subdirectory."))
    }

    positions_file <- NULL
    for (pattern in config@file_patterns$positions) {
        candidate <- file.path(spatial_dir, pattern)
        if (file.exists(candidate)) {
            positions_file <- candidate
            break
        }
    }

    if (is.null(positions_file)) {
        stop(wmsg("No tissue_positions file found in: ", spatial_dir, ". ",
                  "Tried: ", paste(config@file_patterns$positions, collapse = ", ")))
    }

    scalefactors_file <- file.path(spatial_dir,
                                   config@file_patterns$scalefactors)
    if (!file.exists(scalefactors_file)) {
        stop(wmsg("No scalefactors_json.json found in: ", spatial_dir))
    }

    image_files <- list.files(spatial_dir,
                              pattern = config@file_patterns$images,
                              full.names = TRUE)

    list(matrix = matrix_path,
         positions = positions_file,
         scalefactors = scalefactors_file,
         spatial_dir = spatial_dir,
         images = image_files)
}

## Technology-specific helpers retained below

.convert_positions_to_micron <- function(positions, scale, config) {    
    spot_diameter_px <- scale$spot_diameter_fullres
    spot_diameter_um <- config@spot_diameter
    um_per_pixel <- spot_diameter_um / spot_diameter_px

    positions$pxl_row_in_fullres <- positions$pxl_row_in_fullres * um_per_pixel
    positions$pxl_col_in_fullres <- positions$pxl_col_in_fullres * um_per_pixel

    positions
}

.get_spot_radius <- function(scale, unit, config) {
    if (unit == "pixel") {
        radius <- scale$spot_diameter_fullres / 2
    } else {
        radius <- config@spot_diameter / 2
    }
    radius
}

#' @importFrom S4Vectors DataFrame
#' @importFrom MultiAssayExperiment ExperimentList
.build_mase_from_visium <-
function(counts, positions, geometries, imgData, sample_id)
{    
    barcodes <- colnames(counts)

    experiments <- ExperimentList(visium = counts)

    colData <- DataFrame(sample_id = sample_id, row.names = barcodes)
    
    sampleMap <- DataFrame(assay = "visium",
                           primary = barcodes,
                           colname = barcodes)
    
    ## Map positions cell_id to barcodes (component reader standardizes names)
    ## If cell_id column exists from component reader, use it; otherwise try barcode
    barcode_col <- if ("cell_id" %in% colnames(positions)) "cell_id" else "barcode"

    geom_barcodes <- as.character(geometries$instance_id)
    spatialMap <- DataFrame(assay = "visium",
                            colname = barcodes,
                            element_type = "shapes",
                            region = "spots",
                            instance_id = barcodes)
    spatialMap <- spatialMap[spatialMap$instance_id %in% geom_barcodes, , 
                             drop = FALSE]

    shapes <- ShapesLayerList(spots = geometries)

    MultiAssaySpatialExperiment(experiments = experiments,
                                colData = colData,
                                sampleMap = sampleMap,
                                shapes = shapes,
                                imgData = imgData,
                                spatialMap = spatialMap)
}
