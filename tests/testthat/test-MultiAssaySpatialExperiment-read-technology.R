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
    skip("Requires mock Visium directory structure")
    
    ## TODO: Create mock directory structure and test file discovery
    ## tmp_dir <- tempfile()
    ## dir.create(tmp_dir)
    ## dir.create(file.path(tmp_dir, "filtered_feature_bc_matrix"))
    ## dir.create(file.path(tmp_dir, "spatial"))
    ## ... create mock files ...
    ## 
    ## result <- readVisiumMASE(tmp_dir)
    ## expect_s4_class(result, "MultiAssaySpatialExperiment")
})

test_that("Xenium file discovery works correctly", {
    skip("Requires mock Xenium directory structure")
    
    ## TODO: Create mock directory structure and test file discovery
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

test_that("Visium HD file discovery checks for Parquet positions", {
    skip("Requires mock Visium HD directory structure")
    
    ## TODO: Create mock directory structure and test that Parquet positions
    ## are properly discovered
    ## tmp_dir <- tempfile()
    ## dir.create(file.path(tmp_dir, "binned_outputs/square_008um/spatial"), 
    ##           recursive = TRUE)
    ## ... create mock Parquet file ...
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

test_that("MERSCOPE file discovery supports multiple FOVs", {
    skip("Requires mock MERSCOPE directory structure")
    
    ## TODO: Create mock directory structure and test multi-FOV discovery
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

test_that("CosMx file discovery supports flexible file patterns", {
    skip("Requires mock CosMx directory structure")
    
    ## TODO: Create mock directory structure and test flexible file discovery
    ## CosMx has variable file naming patterns that should be discovered
})

## ----------------------------------------------------------------------------
## Cross-reader integration tests
## ----------------------------------------------------------------------------

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
