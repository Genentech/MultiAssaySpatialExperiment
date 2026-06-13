## Shared mock vendor directories for integration tests

.mock_cell_square <- function(cell_id, x0, y0, size = 5) {
    data.frame(
        cell_id = rep(cell_id, 5L),
        vertex_x = c(x0, x0 + size, x0 + size, x0, x0),
        vertex_y = c(y0, y0, y0 + size, y0 + size, y0),
        stringsAsFactors = FALSE
    )
}

.write_10x_counts_dir <- function(mtx_dir, barcodes, features, mat) {
    if (!requireNamespace("DropletUtils", quietly = TRUE))
        stop("Package 'DropletUtils' required for mock 10x matrices")
    dir.create(mtx_dir, recursive = TRUE, showWarnings = FALSE)
    DropletUtils::write10xCounts(
        x = mat,
        path = mtx_dir,
        gene.id = features,
        barcodes = barcodes,
        overwrite = TRUE
    )
}

make_mock_visium_dir <- function(n_spots = 3L) {
    root <- tempfile("mock_visium_")
    mtx_dir <- file.path(root, "filtered_feature_bc_matrix")
    spatial_dir <- file.path(root, "spatial")
    dir.create(spatial_dir, recursive = TRUE, showWarnings = FALSE)

    barcodes <- paste0("AAACCTGAG", seq_len(n_spots), "-1")
    features <- c("GENE1", "GENE2")
    mat <- Matrix::Matrix(
        rep(c(1L, 2L, 0L), length.out = length(features) * n_spots),
        nrow = length(features),
        ncol = n_spots,
        dimnames = list(features, barcodes),
        sparse = TRUE
    )
    .write_10x_counts_dir(mtx_dir, barcodes, features, mat)

    positions <- data.frame(
        barcode = barcodes,
        in_tissue = rep(1L, n_spots),
        array_row = seq_len(n_spots),
        array_col = seq_len(n_spots),
        pxl_row_in_fullres = seq(100, 100 + 10L * (n_spots - 1L), 10L),
        pxl_col_in_fullres = seq(200, 200 + 10L * (n_spots - 1L), 10L),
        stringsAsFactors = FALSE
    )
    write.csv(positions, file.path(spatial_dir, "tissue_positions.csv"),
              row.names = FALSE)

    jsonlite::write_json(
        list(tissue_hires_scalef = 1.0,
             tissue_lowres_scalef = 0.1,
             spot_diameter_fullres = 89.43),
        file.path(spatial_dir, "scalefactors_json.json"),
        auto_unbox = TRUE
    )

    root
}

make_mock_xenium_dir <- function(n_cells = 3L) {
    root <- tempfile("mock_xenium_")

    cell_ids <- paste0("cell", seq_len(n_cells))
    features <- c("GENE1", "GENE2")
    mat <- Matrix::Matrix(
        rep(c(3L, 1L, 0L), length.out = length(features) * n_cells),
        nrow = length(features),
        ncol = n_cells,
        dimnames = list(features, cell_ids),
        sparse = TRUE
    )
    .write_10x_counts_dir(file.path(root, "cell_feature_matrix"), cell_ids, features, mat)

    cells <- data.frame(
        cell_id = cell_ids,
        x_centroid = seq(10, 10 + 5L * (n_cells - 1L), 5L),
        y_centroid = seq(20, 20 + 5L * (n_cells - 1L), 5L),
        stringsAsFactors = FALSE
    )
    write.csv(cells, file.path(root, "cells.csv"), row.names = FALSE)

    boundaries <- do.call(rbind, Map(
        .mock_cell_square,
        cell_ids,
        seq(5, 5 + 10L * (n_cells - 1L), 10L),
        seq(10, 10 + 10L * (n_cells - 1L), 10L)
    ))
    write.csv(boundaries, file.path(root, "cell_boundaries.csv"), row.names = FALSE)

    root
}

.write_parquet_df <- function(df, path) {
    if (!requireNamespace("arrow", quietly = TRUE))
        stop("Package 'arrow' required for mock Parquet fixtures")
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    arrow::write_parquet(df, path)
}

make_mock_merscope_dir <- function(n_cells = 3L) {
    root <- tempfile("mock_merscope_")
    region <- file.path(root, "region_0")
    dir.create(region, recursive = TRUE, showWarnings = FALSE)

    cell_ids <- paste0("cell", seq_len(n_cells))
    meta <- data.frame(
        cell_id = cell_ids,
        center_x = seq(10, 10 + 5L * (n_cells - 1L), 5L),
        center_y = seq(20, 20 + 5L * (n_cells - 1L), 5L),
        stringsAsFactors = FALSE
    )
    write.csv(meta, file.path(region, "cell_metadata.csv"), row.names = FALSE)

    boundaries <- do.call(rbind, Map(
        .mock_cell_square,
        cell_ids,
        seq(5, 5 + 10L * (n_cells - 1L), 10L),
        seq(10, 10 + 10L * (n_cells - 1L), 10L)
    ))
    write.csv(boundaries, file.path(region, "cell_boundaries.csv"), row.names = FALSE)

    root
}

make_mock_cosmx_dir <- function(n_cells = 3L) {
    root <- tempfile("mock_cosmx_")
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
    cell_ids <- paste0("c", seq_len(n_cells))

    meta <- data.frame(
        cell_id = cell_ids,
        x_global_px = seq(10, 10 + 5L * (n_cells - 1L), 5L),
        y_global_px = seq(20, 20 + 5L * (n_cells - 1L), 5L),
        stringsAsFactors = FALSE
    )
    write.csv(meta, file.path(root, "cell_metadata.csv"), row.names = FALSE)

    boundaries <- do.call(rbind, Map(
        .mock_cell_square,
        cell_ids,
        seq(5, 5 + 10L * (n_cells - 1L), 10L),
        seq(10, 10 + 10L * (n_cells - 1L), 10L)
    ))
    write.csv(boundaries, file.path(root, "cell_boundaries.csv"), row.names = FALSE)

    expr <- data.frame(Gene = c("GENE1", "GENE2"), stringsAsFactors = FALSE)
    for (id in cell_ids)
        expr[[id]] <- c(1L, 2L)
    write.csv(expr, file.path(root, "expression_matrix.csv"), row.names = FALSE)

    root
}

make_mock_visiumhd_dir <- function(n_barcodes = 3L, bin_size = "008") {
    root <- tempfile("mock_visiumhd_")
    bin_dir <- file.path(root, paste0("binned_outputs/square_", bin_size, "um"))
    mtx_dir <- file.path(bin_dir, "filtered_feature_bc_matrix")
    spatial_dir <- file.path(bin_dir, "spatial")
    dir.create(spatial_dir, recursive = TRUE, showWarnings = FALSE)

    barcodes <- paste0("BC", seq_len(n_barcodes))
    features <- c("GENE1", "GENE2")
    mat <- Matrix::Matrix(
        rep(c(1L, 2L, 0L), length.out = length(features) * n_barcodes),
        nrow = length(features),
        ncol = n_barcodes,
        dimnames = list(features, barcodes),
        sparse = TRUE
    )
    .write_10x_counts_dir(mtx_dir, barcodes, features, mat)

    positions <- data.frame(
        barcode = barcodes,
        in_tissue = rep(1L, n_barcodes),
        pxl_row_in_fullres = seq(100, 100 + 10L * (n_barcodes - 1L), 10L),
        pxl_col_in_fullres = seq(200, 200 + 10L * (n_barcodes - 1L), 10L),
        stringsAsFactors = FALSE
    )
    .write_parquet_df(positions,
                      file.path(spatial_dir, "tissue_positions.parquet"))

    root
}
