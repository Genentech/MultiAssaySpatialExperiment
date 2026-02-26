### =========================================================================
### Component readers for MultiAssayExperiment slots
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### ColData readers (cell/spot metadata)
###

.clean_coldata_names <- function(df, technology) {
    config <- get_technology_config(technology)
    cols <- find_core_coord_cols(df, config = config)

    if (!is.null(cols$id)) {
        colnames(df)[colnames(df) == cols$id] <- "cell_id"
    }
    if (!is.null(cols$x)) {
        colnames(df)[colnames(df) == cols$x] <- "x_centroid"
    }
    if (!is.null(cols$y)) {
        colnames(df)[colnames(df) == cols$y] <- "y_centroid"
    }

    df
}

.read_coldata_parquet <- function(parquet_path, technology, ...) {
    df <- readParquetForMASE(parquet_path, ...)
    .clean_coldata_names(df, technology)
}

.read_coldata_csv <- function(csv_path, technology, ...) {
    df <- readCSVForMASE(csv_path, ...)
    .clean_coldata_names(df, technology)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Metadata readers
###

#' @importFrom S4Vectors wmsg
.read_visium_scalefactors <- function(json_path) {
    if (!requireNamespace("jsonlite", quietly = TRUE))
        stop(wmsg("Package 'jsonlite' required for scale factors reading. ",
                  "Install with: install.packages('jsonlite')"))

    if (!file.exists(json_path)) {
        warning(wmsg("Scale factors file not found: ", json_path,
                    ". Using defaults."))
        return(list(tissue_hires_scalef = 1.0,
                    tissue_lowres_scalef = 1.0,
                    spot_diameter_fullres = 89.43))
    }

    scale <- jsonlite::read_json(json_path,
                                 simplifyVector = TRUE,
                                 simplifyDataFrame = FALSE,
                                 simplifyMatrix = FALSE)

    expected <- c("tissue_hires_scalef", "tissue_lowres_scalef", 
                  "spot_diameter_fullres")
    missing <- setdiff(expected, names(scale))

    if (length(missing) > 0L) {
        warning(wmsg("Missing scale factor fields: ",
                     paste(missing, collapse = ", ")))
    }

    scale
}

#' @importFrom S4Vectors wmsg
.read_metadata_csv <- function(csv_path, key_col = NULL, value_col = NULL, ...) {
    if (!file.exists(csv_path)) {
        warning(wmsg("Metadata CSV not found: ", csv_path))
        return(NULL)
    }

    df <- readCSVForMASE(csv_path, ...)

    if (!is.null(key_col) && !is.null(value_col)) {
        if (!key_col %in% colnames(df) || !value_col %in% colnames(df)) {
            stop(wmsg("Specified key or value columns not found in CSV"))
        }

        metadata <- as.list(df[[value_col]])
        names(metadata) <- df[[key_col]]
        return(metadata)
    }

    df
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Shapes readers (spatial geometries)
###

## Internal dispatch table for automatic format detection
## Pattern-based dispatch following file naming conventions
.SHAPES_FORMAT_DISPATCH <- list(
    list(pattern = "_sf\\.parquet$", reader = ".read_geometries_geoparquet", priority = 10),
    list(pattern = "\\.parquet$", reader = ".read_geometries_vertex_parquet", priority = 5),
    list(pattern = "\\.geojson$", reader = ".read_geometries_geojson", priority = 5),
    list(pattern = "\\.csv$", reader = ".read_geometries_vertex_csv", priority = 3)
)

## Dispatch to appropriate reader based on file name
.dispatch_shapes_reader <- function(file_path) {
    fname <- basename(file_path)

    matches <- Filter(function(rule) grepl(rule$pattern, fname),
                      .SHAPES_FORMAT_DISPATCH)

    if (length(matches) == 0L) {
        stop(wmsg("No shapes reader found for file: ", file_path))
    }

    ## Sort by priority (highest first) and return top match
    matches <- matches[order(sapply(matches, `[[`, "priority"),
                             decreasing = TRUE)]
    matches[[1L]]$reader
}

.read_geometries_auto <- function(file_path, layer_name = "cells", ...) {
    reader_name <- .dispatch_shapes_reader(file_path)
    reader_fun <- get(reader_name)
    reader_fun(file_path, layer_name = layer_name, ...)
}

.read_geometries_geoparquet <-
function(parquet_path, layer_name = "cells", min_area = NULL, ...) {
    sf_obj <- readGeoParquetForMASE(parquet_path, ...)

    if (!is.null(min_area)) {
        sf_obj <- .filter_polygons(sf_obj, min_area = min_area)
    }

    .check_st_valid(sf_obj, repair = TRUE)
}

.read_geometries_vertex_parquet <-
function(parquet_path, layer_name = "cells",
         cell_id_col = c("cell_id", "Cell"),
         x_col = c("vertex_x", "x_location"),
         y_col = c("vertex_y", "y_location"),
         min_area = NULL, ...)
{    
    df <- readParquetForMASE(parquet_path, ...)
    df <- as.data.frame(df)

    id_actual <- find_column(df, cell_id_col, required = TRUE)
    x_actual <- find_column(df, x_col, required = TRUE)
    y_actual <- find_column(df, y_col, required = TRUE)

    sf_obj <- .vertices_to_polygons(df,
                                    x_col = x_actual,
                                    y_col = y_actual,
                                    id_col = id_actual)

    if (!is.null(min_area)) {
        sf_obj <- .filter_polygons(sf_obj, min_area = min_area)
    }

    .check_st_valid(sf_obj, repair = TRUE)
}

.read_geometries_geojson <-
function(geojson_path, layer_name = "cells", min_area = NULL, ...) {
    sf_obj <- readGeoJSONForMASE(geojson_path, ...)

    if (!is.null(min_area)) {
        sf_obj <- .filter_polygons(sf_obj, min_area = min_area)
    }

    .check_st_valid(sf_obj, repair = TRUE)
}

.read_geometries_vertex_csv <-
function(csv_path, layer_name = "cells",
         cell_id_col = c("cell_id", "Cell"),
         x_col = c("vertex_x", "x_location"),
         y_col = c("vertex_y", "y_location"),
         min_area = NULL, ...)
{
    df <- readCSVForMASE(csv_path, ...)

    id_actual <- find_column(df, cell_id_col, required = TRUE)
    x_actual <- find_column(df, x_col, required = TRUE)
    y_actual <- find_column(df, y_col, required = TRUE)

    sf_obj <- .vertices_to_polygons(df,
                                    x_col = x_actual,
                                    y_col = y_actual,
                                    id_col = id_actual)

    if (!is.null(min_area)) {
        sf_obj <- .filter_polygons(sf_obj, min_area = min_area)
    }

    .check_st_valid(sf_obj, repair = TRUE)
}

## Generates circular polygons from spot centers
#' @importFrom sf st_as_sf st_buffer
.create_spot_geometries <-
function(positions_df, x_col = "x_centroid", y_col = "y_centroid",
         radius = NULL, min_area = NULL)
{
    if (!x_col %in% colnames(positions_df)) {
        stop(wmsg("Column '", x_col, "' not found in positions data"))
    }
    if (!y_col %in% colnames(positions_df)) {
        stop(wmsg("Column '", y_col, "' not found in positions data"))
    }

    centers <- positions_df
    names(centers)[names(centers) == x_col] <- "x"
    names(centers)[names(centers) == y_col] <- "y"

    pts <- st_as_sf(centers, coords = c("x", "y"))
    spots <- st_buffer(pts, dist = radius, nQuadSegs = 16L)

    if (!is.null(min_area)) {
        spots <- .filter_polygons(spots, min_area = min_area)
    }

    spots
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Points readers (transcript coordinates)
###

.read_transcripts_parquet <-
function(parquet_path,
         x_col = c("x_location", "x"),
         y_col = c("y_location", "y"),
         z_col = c("z_location", "z"),
         feature_col = c("feature_name", "gene", "target"),
         cell_col = c("cell_id", "Cell", "CellId"),
         qc_col = c("qv", "quality", "qscore"),
         ...)
{
    df <- readParquetForMASE(parquet_path, ...)

    x_actual <- find_column(df, x_col, required = TRUE)
    y_actual <- find_column(df, y_col, required = TRUE)
    z_actual <- find_column(df, z_col, required = FALSE)
    feature_actual <- find_column(df, feature_col, required = TRUE)
    cell_actual <- find_column(df, cell_col, required = FALSE)
    qc_actual <- find_column(df, qc_col, required = FALSE)

    names(df)[names(df) == x_actual] <- "x"
    names(df)[names(df) == y_actual] <- "y"
    if (!is.null(z_actual)) {
        names(df)[names(df) == z_actual] <- "z"
    }
    names(df)[names(df) == feature_actual] <- "feature"
    if (!is.null(cell_actual)) {
        names(df)[names(df) == cell_actual] <- "cell_id"
    }
    if (!is.null(qc_actual)) {
        names(df)[names(df) == qc_actual] <- "qc_score"
    }

    df
}

.read_transcripts_csv <-
function(csv_path, x_col = c("x", "x_location"), y_col = c("y", "y_location"),
         feature_col = c("gene", "target", "feature_name"), ...)
{
    df <- readCSVForMASE(csv_path, ...)

    x_actual <- find_column(df, x_col, required = TRUE)
    y_actual <- find_column(df, y_col, required = TRUE)
    feature_actual <- find_column(df, feature_col, required = TRUE)

    names(df)[names(df) == x_actual] <- "x"
    names(df)[names(df) == y_actual] <- "y"
    names(df)[names(df) == feature_actual] <- "feature"

    df
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Image metadata readers
###

#' @importFrom S4Vectors wmsg
.read_image_metadata <-
function(image_path, sample_id, image_id, scale_factors = NULL,
         load_image = FALSE)
{
    if (!file.exists(image_path)) {
        warning(wmsg("Image file not found: ", image_path, ". Skipping."))
        return(NULL)
    }

    .get_imgData(image_path = image_path, sample_id = sample_id,
                 image_id = image_id, scale_factors = scale_factors,
                 load_image = load_image)
}

#' @importFrom S4Vectors wmsg
#' @importFrom tools file_path_sans_ext
.read_images_metadata <-
function(image_paths, sample_id, image_ids = NULL, scale_factors = NULL,
         load_images = FALSE)
{
    if (is.null(image_paths) || length(image_paths) == 0L) {
        return(NULL)
    }

    if (is.null(image_ids)) {
        image_ids <- basename(image_paths)
        image_ids <- file_path_sans_ext(image_ids)
    }

    if (length(image_paths) != length(image_ids)) {
        stop(wmsg("Length of image_paths and image_ids must match"))
    }

    imgdata_list <- lapply(seq_along(image_paths), function(i) {
        .read_image_metadata(image_path = image_paths[i],
                             sample_id = sample_id,
                             image_id = image_ids[i],
                             scale_factors = scale_factors,
                             load_image = load_images)
    })

    imgdata_list <- Filter(Negate(is.null), imgdata_list)

    if (length(imgdata_list) == 0L) {
        return(NULL)
    }

    do.call(rbind, imgdata_list)
}

#' @importFrom S4Vectors wmsg
#' @importFrom tools file_path_sans_ext
.read_visium_images <-
function(spatial_dir, sample_id, load_images = FALSE, scale_factors = NULL) {    
    image_files <- c("tissue_lowres_image.png", "tissue_hires_image.png",
                     "aligned_fiducials.jpg")

    image_paths <- file.path(spatial_dir, image_files)

    exists <- file.exists(image_paths)
    if (!any(exists)) {
        warning(wmsg("No image files found in: ", spatial_dir))
        return(NULL)
    }

    image_paths <- image_paths[exists]
    image_ids <- file_path_sans_ext(basename(image_paths))

    .read_images_metadata(image_paths = image_paths,
                          sample_id = sample_id,
                          image_ids = image_ids,
                          scale_factors = scale_factors,
                          load_images = load_images)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Label raster readers
###

#' @importFrom S4Vectors wmsg
.read_label_raster <- function(tiff_path, layer_name = "segmentation", ...) {
    if (!file.exists(tiff_path)) {
        warning(wmsg("Label file not found: ", tiff_path, ". Skipping."))
        return(NULL)
    }

    if (!requireNamespace("terra", quietly = TRUE)) {
        stop(wmsg("Package 'terra' required for label rasters. ",
                  "Install with: BiocManager::install('terra')"))
    }

    r <- terra::rast(tiff_path, ...)
    names(r) <- layer_name

    r
}

#' @importFrom S4Vectors wmsg
.read_label_rasters <- function(tiff_paths, layer_names = NULL, ...) {
    if (is.null(tiff_paths) || length(tiff_paths) == 0L) {
        return(NULL)
    }

    if (is.null(layer_names)) {
        layer_names <- paste0("segmentation_", seq_along(tiff_paths))
    }

    if (length(tiff_paths) != length(layer_names)) {
        stop(wmsg("Length of tiff_paths and layer_names must match"))
    }

    rasters <- lapply(seq_along(tiff_paths), function(i) {
        .read_label_raster(tiff_path = tiff_paths[i],
                           layer_name = layer_names[i],
                           ...)
    })

    rasters <- Filter(Negate(is.null), rasters)

    if (length(rasters) == 0L) {
        return(NULL)
    }

    if (length(rasters) == 1L) {
        rasters[[1L]]
    } else {
        rasters
    }
}

#' @importFrom S4Vectors wmsg
.read_label_ometiff <-
function(ometiff_path, series = 1L, layer_name = "segmentation", ...) {
    if (!file.exists(ometiff_path)) {
        warning(wmsg("OME-TIFF file not found: ", ometiff_path, ". Skipping."))
        return(NULL)
    }

    if (!requireNamespace("RBioFormats", quietly = TRUE)) {
        stop(wmsg("Package 'RBioFormats' required for OME-TIFF. ",
                  "Install with: BiocManager::install('RBioFormats')"))
    }

    img <- RBioFormats::read.image(ometiff_path, series = series, ...)

    if (!inherits(img, "SpatRaster")) {
        if (!requireNamespace("terra", quietly = TRUE)) {
            stop(wmsg("Package 'terra' required for OME-TIFF reading. ",
                      "Install with: BiocManager::install('terra')"))
        }
        img <- terra::rast(img)
    }

    names(img) <- layer_name

    img
}
