### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constants
###

VALID_TECHNOLOGIES <- c("Visium", "VisiumHD", "Xenium", "Vizgen", "CosMx", "MASE")

VALID_DATA_TYPES <- c("filtered", "raw")

VALID_MATRIX_FORMATS <- c("HDF5", "sparse")

VALID_COORD_UNITS <- c("pixel", "micron")

DEFAULT_BLOCK_SIZE <- 100 * 1024 * 1024  # 100 MB (following DelayedArray)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Basic utilities
###

get_block_size <- function() {
    getOption("MASE.block.size", DEFAULT_BLOCK_SIZE)
}

set_block_size <- function(size) {
    if (!isSingleNumber(size) || size <= 0L)
        stop(wmsg("Block size must be a positive number (bytes)"))

    options(MASE.block.size = size)
    invisible(NULL)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### MASE-specific normalization functions
###

#' @importFrom S4Vectors isSingleString wmsg
normarg_data_dir <- function(data_dir) {
    if (!isSingleString(data_dir))
        stop(wmsg("'data_dir' must be a single string"))

    if (!dir.exists(data_dir))
        stop(wmsg("'data_dir' not found: ", data_dir))

    data_dir
}

#' @importFrom S4Vectors isSingleString wmsg
normarg_sample_id <-
function(sample_id, data_dir = NULL, default = "sample_01") {
    if (is.null(sample_id)) {
        if (!is.null(data_dir) && isSingleString(data_dir)) {
            sample_id <- basename(data_dir)
        } else {
            sample_id <- default
        }
    }

    if (!isSingleString(sample_id))
        stop(wmsg("'sample_id' must be NULL or a single string"))

    sample_id
}

#' @importFrom S4Vectors isSingleString wmsg
normarg_technology <-
function(technology, valid_technologies = VALID_TECHNOLOGIES) {
    if (!isSingleString(technology))
        stop(wmsg("'technology' must be a single string"))

    match.arg(technology, valid_technologies)
}

#' @importFrom S4Vectors isSingleString wmsg
normarg_data_type <- function(data, valid_types = VALID_DATA_TYPES) {
    match.arg(data, valid_types)
}

#' @importFrom S4Vectors isSingleString wmsg
normarg_matrix_format <- function(type, valid_formats = VALID_MATRIX_FORMATS) {
    match.arg(type, valid_formats)
}

#' @importFrom S4Vectors isSingleString wmsg
normarg_unit <- function(unit, valid_units = VALID_COORD_UNITS) {
    match.arg(unit, valid_units)
}

.observation_ids <- function(counts) {
    cn <- colnames(counts)
    if (!is.null(cn) && length(cn) > 0L)
        return(cn)
    if (is(counts, "SummarizedExperiment")) {
        cd <- colData(counts)
        if ("Barcode" %in% colnames(cd))
            return(as.character(cd[["Barcode"]]))
        cn <- colnames(cd)
        if (!is.null(cn) && length(cn) > 0L)
            return(cn)
        rn <- rownames(cd)
        if (!is.null(rn) && length(rn) > 0L)
            return(rn)
    }
    character(0L)
}

#' @importFrom S4Vectors isSingleNumber wmsg
normarg_min_area <- function(min_area) {
    if (is.null(min_area))
        return(NULL)

    if (!isSingleNumber(min_area) || min_area < 0L)
        stop(wmsg("'min_area' must be NULL or a non-negative number"))

    as.numeric(min_area)
}

#' @importFrom S4Vectors isSingleNumber wmsg
normarg_min_phred <- function(min_phred) {
    if (is.null(min_phred))
        return(NULL)

    if (!isSingleNumber(min_phred) || min_phred < 0L)
        stop(wmsg("'min_phred' must be NULL or a non-negative number"))

    as.integer(min_phred)
}

#' @importFrom S4Vectors isSingleNumber wmsg
normarg_block_size <- function(block_size) {
    if (is.null(block_size))
        return(get_block_size())

    if (!isSingleNumber(block_size) || block_size <= 0L)
        stop(wmsg("'block_size' must be NULL or a positive number (bytes)"))

    as.numeric(block_size)
}

#' @importFrom S4Vectors isTRUEorFALSE wmsg
normarg_images <- function(images, load_images) {
    if (!isTRUEorFALSE(images))
        stop(wmsg("'images' must be TRUE or FALSE"))

    if (!isTRUEorFALSE(load_images))
        stop(wmsg("'load_images' must be TRUE or FALSE"))

    list(images = images, load_images = load_images)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Column discovery
###

#' @importFrom S4Vectors wmsg
.normarg_field <- function(field, what) {
    if (!is.character(field) || any(is.na(field)))
        stop(wmsg("'", what, "' must be a character vector with no NAs"))
    if (length(field) == 0L)
        stop(wmsg("'", what, "' cannot be empty"))
    field
}

#' @importFrom S4Vectors wmsg
find_column <-
function(df, candidates, required = FALSE, context = "data frame") {
    if (!is.character(candidates) || length(candidates) == 0L)
        stop(wmsg("'candidates' must be a non-empty character vector"))

    df_cols <- colnames(df)
    matches <- candidates %in% df_cols

    if (!any(matches)) {
        if (required) {
            stop(wmsg(
                "Could not find required column in ", context, ". ",
                "Tried: ", paste(candidates, collapse = ", "), ". ",
                "Available columns: ", paste(df_cols, collapse = ", ")
            ))
        }
        return(NULL)
    }

    candidates[matches][1L]
}

find_core_coord_cols <-
function(df, x.field = NULL, y.field = NULL, id.field = NULL, config = NULL) {
    if (!is.null(config)) {
        if (is.null(x.field))
            x.field <- config@column_candidates$x
        if (is.null(y.field))
            y.field <- config@column_candidates$y
        if (is.null(id.field)) {
            id.field <- if (!is.null(config@column_candidates$barcode)) {
                config@column_candidates$barcode
            } else if (!is.null(config@column_candidates$cell_id)) {
                config@column_candidates$cell_id
            } else {
                config@column_candidates[[1L]]
            }
        }
    }

    ## Validate field arguments
    x.field <- .normarg_field(x.field, "x.field")
    y.field <- .normarg_field(y.field, "y.field")
    if (!is.null(id.field))
        id.field <- .normarg_field(id.field, "id.field")

    ## Find columns
    x_col <- find_column(df, x.field, required = TRUE,
                         context = "coordinate data")
    y_col <- find_column(df, y.field, required = TRUE,
                         context = "coordinate data")
    id_col <- if (!is.null(id.field)) {
        find_column(df, id.field, required = FALSE,
                    context = "coordinate data")
    } else {
        NULL
    }

    list(x = x_col, y = y_col, id = id_col)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Geometry utilities
###

#' @importFrom S4Vectors wmsg
#' @importFrom sf st_is_valid st_area
.filter_polygons <- function(geom, min_area = NULL, remove_invalid = TRUE) {
    if (!inherits(geom, "sf") && !inherits(geom, "sfc"))
        stop(wmsg("'geom' must be an sf or sfc object"))

    if (remove_invalid) {
        valid <- st_is_valid(geom)
        if (!all(valid)) {
            message(sprintf("Removing %d invalid geometries", sum(!valid)))
            geom <- geom[valid, ]
        }
    }

    if (!is.null(min_area) && min_area > 0L) {
        areas <- st_area(geom)
        keep <- areas >= min_area
        if (!all(keep)) {
            message(sprintf("Removing %d geometries with area < %g",
                            sum(!keep), min_area))
            geom <- geom[keep, ]
        }
    }

    geom
}

#' @importFrom S4Vectors wmsg
#' @importFrom sf st_is_valid st_buffer
.check_st_valid <- function(geom, repair = TRUE) {

    valid <- st_is_valid(geom)

    if (all(valid))
        return(geom)

    n_invalid <- sum(!valid)
    message(sprintf("Found %d invalid geometries", n_invalid))

    if (repair) {
        message("Attempting repair with st_buffer(x, 0)...")
        geom[!valid] <- st_buffer(geom[!valid], 0)

        still_invalid <- !st_is_valid(geom)
        if (any(still_invalid)) {
            warning(wmsg(sprintf(
                "%d geometries could not be repaired and will be removed",
                sum(still_invalid)
            )))
            geom <- geom[!still_invalid]
        } else {
            message("All geometries repaired successfully")
        }
    } else {
        stop(wmsg(sprintf("%d invalid geometries found. Set repair=TRUE to fix.",
                          n_invalid)))
    }

    geom
}

#' @importFrom S4Vectors wmsg
#' @importFrom sf st_polygon st_sfc st_sf
.vertices_to_polygons <-
function(df, x_col = "vertex_x", y_col = "vertex_y", id_col = "cell_id") {
    required <- c(x_col, y_col, id_col)
    missing <- setdiff(required, colnames(df))
    if (length(missing) > 0L) {
        stop(wmsg("Vertex data missing columns: ",
                  paste(missing, collapse = ", ")))
    }

    xv <- df[[x_col]]
    yv <- df[[y_col]]
    ids <- df[[id_col]]
    unique_ids <- unique(ids)

    ## Build polygons
    poly_list <- lapply(unique_ids, function(id) {
        idx <- ids == id
        poly_coords <- cbind(xv[idx], yv[idx])

        ## Close polygon if not already closed
        if (nrow(poly_coords) > 0L) {
            first_row <- poly_coords[1L, , drop = TRUE]
            last_row <- poly_coords[nrow(poly_coords), , drop = TRUE]
            if (!all(first_row == last_row)) {
                poly_coords <- rbind(poly_coords, first_row)
            }
        }

        st_polygon(list(poly_coords))
    })

    geom_col <- st_sfc(poly_list)
    st_sf(DataFrame(instance_id = as.character(unique_ids)), geometry = geom_col)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Image utilities
###

#' @importFrom S4Vectors wmsg
.get_image_info <- function(image_path, load_data = FALSE) {
    ext <- tolower(tools::file_ext(image_path))

    if (ext %in% c("tif", "tiff")) {
        if (!requireNamespace("terra", quietly = TRUE)) {
            stop(wmsg("Package 'terra' required for TIFF images. ",
                      "Install with: BiocManager::install('terra')"))
        }

        r <- terra::rast(image_path)
        info <- list(width = terra::ncol(r),
                     height = terra::nrow(r),
                     format = "TIFF",
                     data = if (load_data) r else NULL)
    } else if (ext %in% c("png", "jpg", "jpeg")) {
        if (!requireNamespace("terra", quietly = TRUE)) {
            stop(wmsg("Package 'terra' required for image files. ",
                      "Install with: BiocManager::install('terra')"))
        }

        r <- terra::rast(image_path)
        info <- list(width = terra::ncol(r),
                     height = terra::nrow(r),
                     format = toupper(ext),
                     data = if (load_data) r else NULL)
        
    } else if (grepl("ome\\.tif", image_path)) {
        if (!requireNamespace("RBioFormats", quietly = TRUE)) {
            stop(wmsg("Package 'RBioFormats' required for OME-TIFF images. ",
                      "Install with: BiocManager::install('RBioFormats')"))
        }

        img <- RBioFormats::read.image(image_path)
        info <- list(width = dim(img)[1L],
                     height = dim(img)[2L],
                     format = "OME-TIFF",
                     data = if (load_data) img else NULL)
    } else {
        stop(wmsg("Unsupported image format: ", ext, ". ",
                  "Supported: png, jpg, jpeg, tif, tiff, ome.tif"))
    }

    info
}

#' @importFrom S4Vectors DataFrame wmsg
.get_imgData <-
function(image_path,
         sample_id,
         image_id,
         scale_factors = NULL,
         load_image = FALSE)
{
    if (!file.exists(image_path))
        stop(wmsg("Image file not found: ", image_path))

    img_info <- .get_image_info(image_path, load_image)

    payload <- if (load_image) img_info$data else NULL
    DataFrame(sample_id = sample_id,
              image_id = image_id,
              data = I(list(payload)),
              scaleFactor = if (!is.null(scale_factors)) scale_factors else 1.0,
              width = img_info$width,
              height = img_info$height,
              path = image_path)
}
