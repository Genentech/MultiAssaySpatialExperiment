### =========================================================================
### S4 Generics for MASE File Format Readers
### -------------------------------------------------------------------------
###

#' Read Parquet Files for MASE
#'
#' @description
#' Read a Parquet file into a \linkS4class{DataFrame} for use in MASE component
#' assembly or custom readers.
#'
#' @param file_path Character path to Parquet file
#' @param col_select Optional character vector of column names to select
#' @param ... Additional arguments passed to format-specific methods
#'
#' @return DataFrame
#'
#' @examples
#' if (requireNamespace("arrow", quietly = TRUE)) {
#'     tmp <- tempfile(fileext = ".parquet")
#'     df <- data.frame(cell_id = paste0("C", 1:3), x = 1:3, y = 1:3)
#'     arrow::write_parquet(df, tmp)
#'     readParquetForMASE(tmp)
#'     unlink(tmp)
#' }
#'
#' @export
setGeneric("readParquetForMASE",
    function(file_path, ...) standardGeneric("readParquetForMASE"))

#' @rdname readParquetForMASE
#' @importFrom S4Vectors DataFrame wmsg
#' @export
setMethod("readParquetForMASE", "character",
    function(file_path, col_select = NULL, ...) {
        if (!requireNamespace("arrow", quietly = TRUE))
            stop(wmsg("Package 'arrow' required for Parquet reading. ",
                      "Install with: install.packages('arrow')"))

        if (is.null(col_select)) {
            df <- arrow::read_parquet(file_path, ...)
        } else {
            df <- arrow::read_parquet(file_path, col_select = col_select, ...)
        }
        DataFrame(df, check.names = FALSE)
    })

#' Read GeoParquet Files for MASE
#'
#' @description
#' Read a GeoParquet file via \pkg{sfarrow} for geometry layers in MASE readers.
#'
#' @param file_path Character path to GeoParquet file
#' @param ... Additional arguments passed to format-specific methods
#'
#' @return sf object or equivalent object
#'
#' @examples
#' if (requireNamespace("sfarrow", quietly = TRUE) &&
#'     requireNamespace("sf", quietly = TRUE)) {
#'     tmp <- tempfile(fileext = ".parquet")
#'     pt <- sf::st_sfc(sf::st_point(c(1, 2)))
#'     sf_obj <- sf::st_sf(id = 1L, geometry = pt)
#'     sfarrow::st_write_parquet(sf_obj, tmp)
#'     readGeoParquetForMASE(tmp)
#'     unlink(tmp)
#' }
#'
#' @export
setGeneric("readGeoParquetForMASE",
    function(file_path, ...) standardGeneric("readGeoParquetForMASE"))

#' @rdname readGeoParquetForMASE
#' @importFrom S4Vectors wmsg
#' @export
setMethod("readGeoParquetForMASE", "character",
    function(file_path, ...) {
        if (!requireNamespace("sfarrow", quietly = TRUE))
            stop(wmsg("Package 'sfarrow' required for GeoParquet reading. ",
                      "Install with: install.packages('sfarrow')"))

        sfarrow::st_read_parquet(file_path, ...)
    })

#' Read HDF5 Files for MASE
#'
#' @description
#' Read an HDF5 matrix file (10x or AnnData/h5ad) into a
#' \linkS4class{SummarizedExperiment} for MASE assembly.
#'
#' @param file_path Character path to HDF5 file
#' @param type Character: "10x" for 10x HDF5, "anndata" for h5ad
#' @param ... Additional arguments passed to format-specific methods
#'
#' @return SummarizedExperiment
#'
#' @examples
#' \dontrun{
#' ## Requires DropletUtils and a 10x filtered_feature_bc_matrix.h5 file
#' readHDF5ForMASE("filtered_feature_bc_matrix.h5", type = "10x")
#' }
#'
#' @export
setGeneric("readHDF5ForMASE",
    function(file_path, type = c("10x", "anndata"), ...) 
        standardGeneric("readHDF5ForMASE"))

#' @rdname readHDF5ForMASE
#' @export
setMethod("readHDF5ForMASE", c("character", "character"),
    function(file_path, type = c("10x", "anndata"), ...) {
        type <- match.arg(type)
        switch(type,
            "10x" = {
                if (!requireNamespace("DropletUtils", quietly = TRUE))
                    stop(wmsg("Package 'DropletUtils' required for 10x HDF5 reading. ",
                              "Install with: BiocManager::install('DropletUtils')"))
                DropletUtils::read10xCounts(file_path, ...)
            },
            "anndata" = {
                if (!requireNamespace("zellkonverter", quietly = TRUE))
                    stop(wmsg("Package 'zellkonverter' required for anndata HDF5 reading. ",
                              "Install with: BiocManager::install('zellkonverter')"))
                ## Convert to SummarizedExperiment
                sce <- zellkonverter::readH5AD(file_path, ...)
                as(sce, "SummarizedExperiment")
            }
        )
    })

#' Read GeoJSON Files for MASE
#'
#' @description
#' Read a GeoJSON file via \pkg{sf} for shape layers in MASE readers.
#'
#' @param file_path Character path to GeoJSON file
#' @param quiet Logical, suppress sf reading messages
#' @param ... Additional arguments passed to format-specific methods
#'
#' @return sf or equivalent object
#'
#' @examples
#' tmp <- tempfile(fileext = ".geojson")
#' pt <- sf::st_point(c(10, 20))
#' sf_obj <- sf::st_sf(id = "C1", geometry = sf::st_sfc(pt))
#' sf::st_write(sf_obj, tmp, quiet = TRUE, delete_dsn = TRUE)
#' readGeoJSONForMASE(tmp, quiet = TRUE)
#' unlink(tmp)
#'
#' @export
setGeneric("readGeoJSONForMASE",
    function(file_path, ...) standardGeneric("readGeoJSONForMASE"))

#' @rdname readGeoJSONForMASE
#' @importFrom sf st_read
#' @export
setMethod("readGeoJSONForMASE", "character",
    function(file_path, quiet = TRUE, ...) {
        st_read(file_path, quiet = quiet, ...)
    })

#' Read CSV Files for MASE
#'
#' @description
#' Read a CSV file into a \linkS4class{DataFrame} for MASE component assembly.
#'
#' @param file_path Character path to CSV file
#' @param ... Additional arguments passed to format-specific methods
#'
#' @return DataFrame
#'
#' @examples
#' tmp <- tempfile(fileext = ".csv")
#' df <- data.frame(
#'     cell_id = paste0("C", 1:3),
#'     x_centroid = 1:3,
#'     y_centroid = 1:3)
#' write.csv(df, tmp, row.names = FALSE)
#' readCSVForMASE(tmp)
#' unlink(tmp)
#'
#' @export
setGeneric("readCSVForMASE",
    function(file_path, ...) standardGeneric("readCSVForMASE"))

#' @rdname readCSVForMASE
#' @importFrom S4Vectors DataFrame
#' @importFrom utils read.csv
#' @export
setMethod("readCSVForMASE", "character",
    function(file_path, ...) {
        DataFrame(read.csv(file_path, ...), check.names = FALSE)
    })
