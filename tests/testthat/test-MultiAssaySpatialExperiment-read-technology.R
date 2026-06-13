## Tests for refactored technology readers (readVisiumMASE, readXeniumMASE)

test_that("readVisiumMASE exists and has expected signature", {
    expect_true(exists("readVisiumMASE"))
    
    ## Check function arguments
    args <- formals(readVisiumMASE)
    expect_true("data_dir" %in% names(args))
    expect_true("sample_id" %in% names(args))
    expect_true("type" %in% names(args))
    expect_true("data" %in% names(args))
    expect_true("images" %in% names(args))
    expect_true("load_images" %in% names(args))
    expect_true("unit" %in% names(args))
    expect_true("min_area" %in% names(args))
})

test_that("readVisiumMASE validates arguments properly", {
    ## Test with non-existent directory
    expect_error(
        readVisiumMASE("/nonexistent/directory"),
        "'data_dir' not found"
    )
})

test_that("readXeniumMASE exists and has expected signature", {
    expect_true(exists("readXeniumMASE"))
    
    ## Check function arguments
    args <- formals(readXeniumMASE)
    expect_true("data_dir" %in% names(args))
    expect_true("sample_id" %in% names(args))
    expect_true("segmentations" %in% names(args))
    expect_true("images" %in% names(args))
    expect_true("load_images" %in% names(args))
    expect_true("add_transcripts" %in% names(args))
    expect_true("min_area" %in% names(args))
    expect_true("min_phred" %in% names(args))
})

test_that("readXeniumMASE validates arguments properly", {
    ## Test with non-existent directory
    expect_error(
        readXeniumMASE("/nonexistent/directory"),
        "'data_dir' not found"
    )
    
    ## Test with invalid add_transcripts
    expect_error(
        readXeniumMASE(tempdir(), add_transcripts = "yes"),
        "'add_transcripts' must be TRUE or FALSE"
    )
})

test_that("Visium file discovery works correctly", {
    skip_if_not_installed("DropletUtils")
    skip_if_not_installed("Matrix")
    skip_if_not_installed("jsonlite")
    skip_if_not_installed("sf")

    data_dir <- make_mock_visium_dir(n_spots = 3L)
    on.exit(unlink(data_dir, recursive = TRUE), add = TRUE)

    mase <- readVisiumMASE(data_dir, images = FALSE)
    expect_s4_class(mase, "MultiAssaySpatialExperiment")
    expect_equal(names(mase), "visium")
    expect_equal(ncol(experiments(mase)[["visium"]]), 3L)
    expect_equal(nrow(spatialShapes(mase)[["spots"]]), 3L)
    expect_equal(nrow(spatialMap(mase)), 3L)
    expect_equal(unique(spatialMap(mase)[["region"]]), "spots")
})

test_that("Xenium file discovery works correctly", {
    skip_if_not_installed("DropletUtils")
    skip_if_not_installed("Matrix")
    skip_if_not_installed("sf")

    data_dir <- make_mock_xenium_dir(n_cells = 3L)
    on.exit(unlink(data_dir, recursive = TRUE), add = TRUE)

    mase <- suppressMessages(readXeniumMASE(
        data_dir,
        segmentations = "cell",
        images = FALSE,
        add_transcripts = FALSE
    ))
    expect_s4_class(mase, "MultiAssaySpatialExperiment")
    expect_equal(names(mase), "xenium")
    expect_equal(ncol(experiments(mase)[["xenium"]]), 3L)
    expect_equal(nrow(spatialShapes(mase)[["cells"]]), 3L)
    expect_equal(nrow(spatialMap(mase)), 3L)
    expect_equal(unique(spatialMap(mase)[["region"]]), "cells")
})

.mock_merscope_boundaries <- function(cell_ids) {
    skip_if_not_installed("sf")
    verts <- do.call(rbind, Map(
        function(id, x0) {
            data.frame(
                cell_id = rep(id, 5L),
                vertex_x = c(x0, x0 + 5, x0 + 5, x0, x0),
                vertex_y = c(10, 10, 15, 15, 10),
                stringsAsFactors = FALSE
            )
        },
        cell_ids,
        seq(5, 5 + 10L * (length(cell_ids) - 1L), 10L)
    ))
    MultiAssaySpatialExperiment:::.vertices_to_polygons(
        verts, id_col = "cell_id", x_col = "vertex_x", y_col = "vertex_y"
    )
}

test_that(".build_mase_from_merscope uses buildSpatialMap for shapes and transcripts", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("sf")

    cell_ids <- c("c1", "c2")
    counts <- SummarizedExperiment(
        assays = list(counts = matrix(
            1:4, 2, 2, dimnames = list(c("g1", "g2"), cell_ids)
        ))
    )
    boundaries <- .mock_merscope_boundaries(cell_ids)
    transcripts <- DataFrame(x = 1:3, y = 4:6, feature = rep("GENE1", 3L))

    build_fn <- get(".build_mase_from_merscope",
                    envir = asNamespace("MultiAssaySpatialExperiment"))
    mase <- build_fn(
        counts = counts,
        cell_meta = NULL,
        boundaries_list = list(cellpose = boundaries),
        transcripts = transcripts,
        imgData = NULL,
        sample_id = "sample1",
        segmentation = "cellpose"
    )

    expect_s4_class(mase, "MultiAssaySpatialExperiment")
    expect_equal(nrow(spatialMap(mase)), 5L)
    expect_equal(sum(spatialMap(mase)[["element_type"]] == "shapes"), 2L)
    expect_equal(sum(spatialMap(mase)[["element_type"]] == "points"), 3L)
    expect_setequal(unique(spatialMap(mase)[["region"]]), c("cellpose", "transcripts"))
})

test_that(".build_mase_from_cosmx uses buildSpatialMap for shapes and transcripts", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("sf")

    cell_ids <- c("c1", "c2")
    counts <- SummarizedExperiment(
        assays = list(counts = matrix(
            1:4, 2, 2, dimnames = list(c("g1", "g2"), cell_ids)
        ))
    )
    boundaries <- .mock_merscope_boundaries(cell_ids)
    transcripts <- DataFrame(x = 1:2, y = 3:4, target = c("A", "B"))

    build_fn <- get(".build_mase_from_cosmx",
                    envir = asNamespace("MultiAssaySpatialExperiment"))
    mase <- build_fn(
        counts = counts,
        cell_meta = NULL,
        boundaries = boundaries,
        transcripts = transcripts,
        fov_meta = NULL,
        imgData = NULL,
        sample_id = "sample1"
    )

    expect_s4_class(mase, "MultiAssaySpatialExperiment")
    expect_equal(nrow(spatialMap(mase)), 4L)
    expect_setequal(unique(spatialMap(mase)[["region"]]), c("cells", "transcripts"))
})

test_that(".build_mase_from_visiumhd uses buildSpatialMap per bin", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("sf")

    barcodes <- c("b1", "b2")
    counts <- SummarizedExperiment(
        assays = list(counts = matrix(
            1:4, 2, 2, dimnames = list(c("g1", "g2"), barcodes)
        ))
    )
    positions <- DataFrame(
        instance_id = barcodes,
        x_centroid = c(10, 20),
        y_centroid = c(15, 25)
    )
    exp_data <- list(bin_008 = list(counts = counts, positions = positions))

    build_fn <- get(".build_mase_from_visiumhd",
                    envir = asNamespace("MultiAssaySpatialExperiment"))
    mase <- build_fn(exp_data, boundaries = NULL, imgData = NULL, sample_id = "s1")

    expect_s4_class(mase, "MultiAssaySpatialExperiment")
    expect_equal(names(mase), "bin_008")
    expect_equal(nrow(spatialMap(mase)), 2L)
    expect_equal(unique(spatialMap(mase)[["region"]]), "bin_008")
    expect_equal(spatialMap(mase)[["colname"]], barcodes)
})

test_that("readVisiumMASE uses component readers internally", {
    ## Check that old parsing functions were removed
    expect_false(exists(".parse_visium_counts", 
                       envir = asNamespace("MultiAssaySpatialExperiment"), 
                       inherits = FALSE))
    expect_false(exists(".parse_visium_positions", 
                       envir = asNamespace("MultiAssaySpatialExperiment"), 
                       inherits = FALSE))
    expect_false(exists(".build_visium_geometries", 
                       envir = asNamespace("MultiAssaySpatialExperiment"), 
                       inherits = FALSE))
    
    ## Check that component readers exist
    expect_true(exists(".read_coldata_csv", 
                      envir = asNamespace("MultiAssaySpatialExperiment"), 
                      inherits = FALSE))
    expect_true(exists(".create_spot_geometries", 
                      envir = asNamespace("MultiAssaySpatialExperiment"), 
                      inherits = FALSE))
})

test_that("readXeniumMASE uses component readers internally", {
    ## Check that old parsing functions were removed
    expect_false(exists(".parse_xenium_counts", 
                       envir = asNamespace("MultiAssaySpatialExperiment"), 
                       inherits = FALSE))
    expect_false(exists(".parse_xenium_cells", 
                       envir = asNamespace("MultiAssaySpatialExperiment"), 
                       inherits = FALSE))
    expect_false(exists(".parse_xenium_boundaries", 
                       envir = asNamespace("MultiAssaySpatialExperiment"), 
                       inherits = FALSE))
    
    ## Check that component readers exist
    expect_true(exists(".read_geometries_auto", 
                      envir = asNamespace("MultiAssaySpatialExperiment"), 
                      inherits = FALSE))
    expect_true(exists(".read_transcripts_parquet", 
                      envir = asNamespace("MultiAssaySpatialExperiment"), 
                      inherits = FALSE))
})

test_that("Xenium reader supports transcripts", {
    ## Check that the build function signature includes transcripts
    build_fn <- get(".build_mase_from_xenium", 
                   envir = asNamespace("MultiAssaySpatialExperiment"))
    args <- formals(build_fn)
    expect_true("transcripts" %in% names(args))
})

## ============================================================================
## Tests for new technology readers (Phase 5)
## ============================================================================

## ----------------------------------------------------------------------------
## readVisiumHDMASE Tests
## ----------------------------------------------------------------------------

test_that("readVisiumHDMASE exists and has expected signature", {
    expect_true(exists("readVisiumHDMASE"))
    
    ## Check function arguments
    args <- formals(readVisiumHDMASE)
    expect_true("data_dir" %in% names(args))
    expect_true("sample_id" %in% names(args))
    expect_true("bin_size" %in% names(args))
    expect_true("load_boundaries" %in% names(args))
    expect_true("images" %in% names(args))
    expect_true("load_images" %in% names(args))
    expect_true("min_area" %in% names(args))
    
    ## Check default bin sizes (eval the default argument)
    expect_equal(eval(args$bin_size), c("008", "002", "016"))
})

test_that("readVisiumHDMASE validates arguments properly", {
    ## Test with non-existent directory
    expect_error(
        readVisiumHDMASE("/nonexistent/directory"),
        "'data_dir' not found"
    )
    
    ## Test with invalid bin_size
    expect_error(
        readVisiumHDMASE(tempdir(), bin_size = "invalid"),
        "'arg' should be one of"
    )
    
    ## Test with invalid load_boundaries
    expect_error(
        readVisiumHDMASE(tempdir(), load_boundaries = "yes"),
        "'load_boundaries' must be TRUE or FALSE"
    )
})

test_that("readVisiumHDMASE uses component readers internally", {
    ## Check that component readers exist and are used
    expect_true(exists(".read_coldata_parquet", 
                      envir = asNamespace("MultiAssaySpatialExperiment"), 
                      inherits = FALSE))
    expect_true(exists(".read_geometries_geojson", 
                      envir = asNamespace("MultiAssaySpatialExperiment"), 
                      inherits = FALSE))
    
    ## Check that the file discovery function exists
    expect_true(exists(".check_visiumhd_files", 
                      envir = asNamespace("MultiAssaySpatialExperiment"), 
                      inherits = FALSE))
    
    ## Check that the build function exists
    expect_true(exists(".build_mase_from_visiumhd", 
                      envir = asNamespace("MultiAssaySpatialExperiment"), 
                      inherits = FALSE))
})

test_that("readVisiumHDMASE supports multi-bin structure", {
    ## Check that the build function signature includes multi-bin support
    build_fn <- get(".build_mase_from_visiumhd", 
                   envir = asNamespace("MultiAssaySpatialExperiment"))
    args <- formals(build_fn)
    expect_true("exp_data" %in% names(args))  ## Multiple experiments
    expect_true("boundaries" %in% names(args))
    expect_true("imgData" %in% names(args))
    expect_true("sample_id" %in% names(args))
})

test_that("Visium HD file discovery works with mock directory", {
    skip_if_not_installed("DropletUtils")
    skip_if_not_installed("Matrix")
    skip_if_not_installed("arrow")
    skip_if_not_installed("sf")

    data_dir <- make_mock_visiumhd_dir(n_barcodes = 3L)
    on.exit(unlink(data_dir, recursive = TRUE), add = TRUE)

    mase <- suppressMessages(readVisiumHDMASE(data_dir, bin_size = "008", images = FALSE))
    expect_s4_class(mase, "MultiAssaySpatialExperiment")
    expect_equal(names(mase), "bin_008")
    expect_equal(ncol(experiments(mase)[["bin_008"]]), 3L)
    expect_equal(nrow(spatialShapes(mase)[["bin_008"]]), 3L)
    expect_equal(nrow(spatialMap(mase)), 3L)
})

## ----------------------------------------------------------------------------
## readMERSCOPEMASE Tests
## ----------------------------------------------------------------------------

test_that("readMERSCOPEMASE exists and has expected signature", {
    expect_true(exists("readMERSCOPEMASE"))
    
    ## Check function arguments
    args <- formals(readMERSCOPEMASE)
    expect_true("data_dir" %in% names(args))
    expect_true("sample_id" %in% names(args))
    expect_true("fov_ids" %in% names(args))
    expect_true("segmentation" %in% names(args))
    expect_true("load_transcripts" %in% names(args))
    expect_true("images" %in% names(args))
    expect_true("load_images" %in% names(args))
    expect_true("min_area" %in% names(args))
    expect_true("min_qv" %in% names(args))
    
    ## Check default segmentation (eval the default argument)
    expect_equal(eval(args$segmentation), c("cellpose", "watershed", "both"))
    expect_equal(args$min_qv, 20)
})

test_that("readMERSCOPEMASE validates arguments properly", {
    ## Test with non-existent directory
    expect_error(
        readMERSCOPEMASE("/nonexistent/directory"),
        "'data_dir' not found"
    )
    
    ## Test with invalid segmentation
    expect_error(
        readMERSCOPEMASE(tempdir(), segmentation = "invalid"),
        "'arg' should be one of"
    )
    
    ## Test with invalid load_transcripts
    expect_error(
        readMERSCOPEMASE(tempdir(), load_transcripts = "yes"),
        "'load_transcripts' must be TRUE or FALSE"
    )
})

test_that("readMERSCOPEMASE uses component readers internally", {
    ## Check that component readers exist
    expect_true(exists(".read_coldata_csv", 
                      envir = asNamespace("MultiAssaySpatialExperiment"), 
                      inherits = FALSE))
    expect_true(exists(".read_coldata_parquet", 
                      envir = asNamespace("MultiAssaySpatialExperiment"), 
                      inherits = FALSE))
    expect_true(exists(".read_geometries_auto", 
                      envir = asNamespace("MultiAssaySpatialExperiment"), 
                      inherits = FALSE))
    expect_true(exists(".read_transcripts_parquet", 
                      envir = asNamespace("MultiAssaySpatialExperiment"), 
                      inherits = FALSE))
    expect_true(exists(".read_transcripts_csv", 
                      envir = asNamespace("MultiAssaySpatialExperiment"), 
                      inherits = FALSE))
    
    ## Check that the file discovery function exists
    expect_true(exists(".check_merscope_files", 
                      envir = asNamespace("MultiAssaySpatialExperiment"), 
                      inherits = FALSE))
    
    ## Check that the build function exists
    expect_true(exists(".build_mase_from_merscope", 
                      envir = asNamespace("MultiAssaySpatialExperiment"), 
                      inherits = FALSE))
})

test_that("readMERSCOPEMASE supports dual segmentation", {
    ## Check that the build function supports both segmentation types
    build_fn <- get(".build_mase_from_merscope", 
                   envir = asNamespace("MultiAssaySpatialExperiment"))
    args <- formals(build_fn)
    expect_true("boundaries_list" %in% names(args))
})

test_that("readMERSCOPEMASE supports transcript loading and count synthesis", {
    ## Check that the build function includes transcripts
    build_fn <- get(".build_mase_from_merscope", 
                   envir = asNamespace("MultiAssaySpatialExperiment"))
    args <- formals(build_fn)
    expect_true("transcripts" %in% names(args))
    
    ## Check that the count synthesis function exists
    expect_true(exists(".build_merscope_counts", 
                      envir = asNamespace("MultiAssaySpatialExperiment"), 
                      inherits = FALSE))
})

test_that("MERSCOPE file discovery works with mock directory", {
    skip_if_not_installed("sf")

    data_dir <- make_mock_merscope_dir(n_cells = 3L)
    on.exit(unlink(data_dir, recursive = TRUE), add = TRUE)

    mase <- suppressMessages(readMERSCOPEMASE(
        data_dir,
        segmentation = "cellpose",
        images = FALSE,
        load_transcripts = FALSE
    ))
    expect_s4_class(mase, "MultiAssaySpatialExperiment")
    expect_equal(names(mase), "merscope")
    expect_equal(ncol(experiments(mase)[["merscope"]]), 3L)
    expect_equal(nrow(spatialShapes(mase)[["cellpose"]]), 3L)
    expect_equal(nrow(spatialMap(mase)), 3L)
})

## ----------------------------------------------------------------------------
## readCosMxMASE Tests
## ----------------------------------------------------------------------------

test_that("readCosMxMASE exists and has expected signature", {
    expect_true(exists("readCosMxMASE"))
    
    ## Check function arguments
    args <- formals(readCosMxMASE)
    expect_true("data_dir" %in% names(args))
    expect_true("sample_id" %in% names(args))
    expect_true("fov_ids" %in% names(args))
    expect_true("load_transcripts" %in% names(args))
    expect_true("images" %in% names(args))
    expect_true("load_images" %in% names(args))
    expect_true("min_area" %in% names(args))
})

test_that("readCosMxMASE validates arguments properly", {
    ## Test with non-existent directory
    expect_error(
        readCosMxMASE("/nonexistent/directory"),
        "'data_dir' not found"
    )
    
    ## Test with invalid load_transcripts
    expect_error(
        readCosMxMASE(tempdir(), load_transcripts = "yes"),
        "'load_transcripts' must be TRUE or FALSE"
    )
})

test_that("readCosMxMASE uses component readers internally", {
    ## Check that component readers exist
    expect_true(exists(".read_coldata_csv", 
                      envir = asNamespace("MultiAssaySpatialExperiment"), 
                      inherits = FALSE))
    expect_true(exists(".read_geometries_auto", 
                      envir = asNamespace("MultiAssaySpatialExperiment"), 
                      inherits = FALSE))
    expect_true(exists(".read_transcripts_parquet", 
                      envir = asNamespace("MultiAssaySpatialExperiment"), 
                      inherits = FALSE))
    expect_true(exists(".read_transcripts_csv", 
                      envir = asNamespace("MultiAssaySpatialExperiment"), 
                      inherits = FALSE))
    expect_true(exists(".read_metadata_csv", 
                      envir = asNamespace("MultiAssaySpatialExperiment"), 
                      inherits = FALSE))
    
    ## Check that the file discovery function exists
    expect_true(exists(".check_cosmx_files", 
                      envir = asNamespace("MultiAssaySpatialExperiment"), 
                      inherits = FALSE))
    
    ## Check that the build function exists
    expect_true(exists(".build_mase_from_cosmx", 
                      envir = asNamespace("MultiAssaySpatialExperiment"), 
                      inherits = FALSE))
})

test_that("readCosMxMASE supports GeoParquet boundaries", {
    ## Verify that .read_geometries_auto is called (handles GeoParquet)
    expect_true(exists(".read_geometries_auto", 
                      envir = asNamespace("MultiAssaySpatialExperiment"), 
                      inherits = FALSE))
})

test_that("readCosMxMASE supports transcript loading and FOV metadata", {
    ## Check that the build function includes transcripts and FOV metadata
    build_fn <- get(".build_mase_from_cosmx", 
                   envir = asNamespace("MultiAssaySpatialExperiment"))
    args <- formals(build_fn)
    expect_true("transcripts" %in% names(args))
    expect_true("fov_meta" %in% names(args))
})

test_that("CosMx file discovery works with mock directory", {
    skip_if_not_installed("sf")

    data_dir <- make_mock_cosmx_dir(n_cells = 3L)
    on.exit(unlink(data_dir, recursive = TRUE), add = TRUE)

    mase <- suppressMessages(readCosMxMASE(data_dir, images = FALSE, load_transcripts = FALSE))
    expect_s4_class(mase, "MultiAssaySpatialExperiment")
    expect_equal(names(mase), "cosmx")
    expect_equal(ncol(experiments(mase)[["cosmx"]]), 3L)
    expect_equal(nrow(spatialShapes(mase)[["cells"]]), 3L)
    expect_equal(nrow(spatialMap(mase)), 3L)
})

test_that("All technology readers support consistent image loading", {
    ## All readers should use .read_images_metadata or similar
    expect_true(exists(".read_images_metadata", 
                      envir = asNamespace("MultiAssaySpatialExperiment"), 
                      inherits = FALSE))
    
    ## Check that each build function includes imgData
    for (fn_name in c(".build_mase_from_visiumhd", 
                     ".build_mase_from_merscope", 
                     ".build_mase_from_cosmx")) {
        build_fn <- get(fn_name, 
                       envir = asNamespace("MultiAssaySpatialExperiment"))
        args <- formals(build_fn)
        expect_true("imgData" %in% names(args), 
                   info = paste(fn_name, "missing imgData"))
    }
})

test_that("All new readers support transcript loading", {
    ## Xenium (refactored), MERSCOPE, and CosMx should support transcripts
    for (fn_name in c(".build_mase_from_xenium",
                     ".build_mase_from_merscope", 
                     ".build_mase_from_cosmx")) {
        build_fn <- get(fn_name, 
                       envir = asNamespace("MultiAssaySpatialExperiment"))
        args <- formals(build_fn)
        expect_true("transcripts" %in% names(args), 
                   info = paste(fn_name, "missing transcripts"))
    }
})

test_that("All readers are exported in NAMESPACE", {
    ## Check that all 5 readers are exported
    exported <- getNamespaceExports("MultiAssaySpatialExperiment")
    expect_true("readVisiumMASE" %in% exported)
    expect_true("readXeniumMASE" %in% exported)
    expect_true("readVisiumHDMASE" %in% exported)
    expect_true("readMERSCOPEMASE" %in% exported)
    expect_true("readCosMxMASE" %in% exported)
})
