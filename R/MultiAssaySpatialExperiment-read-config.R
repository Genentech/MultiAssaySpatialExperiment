### =========================================================================
### Technology configuration registry
### Extensible design for adding new technologies
### -------------------------------------------------------------------------

## Technology configuration class
.TechnologyConfig <- setClass(
    "TechnologyConfig",
    slots = c(
        technology = "character",
        file_patterns = "list",
        column_candidates = "list",
        default_unit = "character",
        spot_diameter = "numeric"
    )
)

## Registry (environment for fast lookup)
.TECHNOLOGY_REGISTRY <- new.env(parent = emptyenv())

## Register a technology
#' @keywords internal
register_technology <- function(config) {
    if (!is(config, "TechnologyConfig"))
        stop(wmsg("'config' must be a TechnologyConfig object"))

    tech <- config@technology
    .TECHNOLOGY_REGISTRY[[tech]] <- config
    invisible(NULL)
}

## Get technology config
#' @keywords internal
get_technology_config <- function(technology) {
    if (!exists(technology, envir = .TECHNOLOGY_REGISTRY, inherits = FALSE))
        stop(wmsg("Technology '", technology, "' not registered"))

    .TECHNOLOGY_REGISTRY[[technology]]
}

## List registered technologies
#' @keywords internal
list_technologies <- function() {
    names(.TECHNOLOGY_REGISTRY)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Pre-defined technology configs
###

## Visium v1/v2
.VISIUM_CONFIG <- new("TechnologyConfig",
    technology = "Visium",
    file_patterns = list(
        matrix = c("filtered_feature_bc_matrix", "raw_feature_bc_matrix"),
        positions = c("tissue_positions.csv", "tissue_positions_list.csv"),
        scalefactors = "scalefactors_json.json",
        images = "\\.(png|tif|tiff)$"
    ),
    column_candidates = list(
        barcode = c("barcode", "barcodes", "cell_id", "spot_id"),
        x = c("pxl_col_in_fullres", "col", "x", "x_coord"),
        y = c("pxl_row_in_fullres", "row", "y", "y_coord"),
        in_tissue = c("in_tissue", "tissue"),
        array_row = c("array_row", "row_index"),
        array_col = c("array_col", "col_index", "column_index")
    ),
    default_unit = "pixel",
    spot_diameter = 55  # microns
)

## Xenium
.XENIUM_CONFIG <- new("TechnologyConfig",
    technology = "Xenium",
    file_patterns = list(
        matrix = c("cell_feature_matrix.h5", "cell_feature_matrix"),
        cells = c("cells.parquet", "cells.csv"),
        cell_boundaries = c("cell_boundaries_sf.parquet", 
                           "cell_boundaries.parquet",
                           "cell_boundaries.csv"),
        nucleus_boundaries = c("nucleus_boundaries_sf.parquet",
                              "nucleus_boundaries.parquet",
                              "nucleus_boundaries.csv"),
        images = "morphology.*\\.ome\\.tif",
        transcripts = c("transcripts.parquet", "detected_transcripts.csv")
    ),
    column_candidates = list(
        cell_id = c("cell_id", "id", "barcode", "label_id"),
        x = c("x_centroid", "x", "x_coord", "center_x"),
        y = c("y_centroid", "y", "y_coord", "center_y"),
        vertex_x = c("vertex_x", "x_location", "x"),
        vertex_y = c("vertex_y", "y_location", "y")
    ),
    default_unit = "micron",
    spot_diameter = NA_real_  # Not applicable
)

## Visium HD
.VISIUMHD_CONFIG <- new("TechnologyConfig",
    technology = "VisiumHD",
    file_patterns = list(
        matrix = c("filtered_feature_bc_matrix.h5", "raw_feature_bc_matrix.h5"),
        positions = "tissue_positions.parquet",  # Parquet format in HD!
        boundaries = c("tissue_positions.geojson", "cell_segmentations.geojson"),
        scalefactors = "scalefactors_json.json",
        images = "\\.(png|tif|tiff)$",
        barcode_mappings = "barcode_mappings.parquet"
    ),
    column_candidates = list(
        barcode = c("barcode", "barcodes", "cell_id"),
        x = c("pxl_col_in_fullres", "x", "x_coord"),
        y = c("pxl_row_in_fullres", "y", "y_coord"),
        in_tissue = c("in_tissue", "tissue")
    ),
    default_unit = "pixel",
    spot_diameter = c(2, 8, 16)  # Multiple bin sizes in microns
)

## MERSCOPE (Vizgen)
.MERSCOPE_CONFIG <- new("TechnologyConfig",
    technology = "MERSCOPE",
    file_patterns = list(
        cells = c("cell_metadata.csv", "cell_metadata.parquet", 
                 "cells.csv", "cells.parquet"),
        cell_boundaries_cellpose = c("cellpose_micron_space.parquet",
                                    "cell_boundaries_cellpose.parquet",
                                    "cell_boundaries.parquet"),
        cell_boundaries_watershed = c("watershed_micron_space.parquet",
                                     "cell_boundaries_watershed.parquet"),
        transcripts = c("detected_transcripts.csv", "transcripts.parquet",
                       "detected_transcripts.parquet"),
        images = "mosaic_.*\\.(tif|tiff|png)$"
    ),
    column_candidates = list(
        cell_id = c("cell_id", "EntityID", "ID"),
        x = c("center_x", "CenterX", "x_global", "x"),
        y = c("center_y", "CenterY", "y_global", "y"),
        vertex_x = c("vertex_x", "x_location", "x"),
        vertex_y = c("vertex_y", "y_location", "y"),
        gene = c("gene", "feature_name", "target"),
        qv = c("qv", "quality", "qscore")
    ),
    default_unit = "micron",
    spot_diameter = NA_real_  # Cell boundaries, not spots
)

## CosMx (NanoString)
.COSMX_CONFIG <- new("TechnologyConfig",
    technology = "CosMx",
    file_patterns = list(
        cells = c("*_metadata_file.csv", "*_cell_metadata.csv", 
                 "cell_metadata.csv"),
        cell_boundaries = c("*_cell_boundaries_sf.parquet",
                           "*_boundaries.parquet",
                           "cell_boundaries.parquet"),
        expression = c("*_exprMat_file.csv", "*_expression.csv",
                      "expression_matrix.csv"),
        transcripts = c("*_tx_file.csv", "*_transcripts.csv",
                       "transcripts.parquet"),
        fov_positions = c("*_fov_positions_file.csv", "fov_positions.csv"),
        images = "\\.(tif|tiff|jpg|jpeg)$"
    ),
    column_candidates = list(
        cell_id = c("cell_ID", "cell_id", "Cell"),
        x = c("CenterX_global_px", "x_global", "x"),
        y = c("CenterY_global_px", "y_global", "y"),
        vertex_x = c("vertex_x", "x_location", "x"),
        vertex_y = c("vertex_y", "y_location", "y"),
        gene = c("target", "gene", "feature_name")
    ),
    default_unit = "pixel",  # Global pixel coordinates
    spot_diameter = NA_real_  # Cell boundaries, not spots
)

## Register configurations on package load
.onLoad_register_technologies <- function() {
    register_technology(.VISIUM_CONFIG)
    register_technology(.XENIUM_CONFIG)
    register_technology(.VISIUMHD_CONFIG)
    register_technology(.MERSCOPE_CONFIG)
    register_technology(.COSMX_CONFIG)
}
