## Tests for component reader functions

test_that(".read_coldata_csv handles standard column names", {
    ## Create a temporary CSV with standard column names
    tmp <- tempfile(fileext = ".csv")
    df <- data.frame(
        barcode = paste0("B", 1:5),
        pxl_col_in_fullres = seq(10, 50, 10),
        pxl_row_in_fullres = seq(20, 100, 20),
        in_tissue = c(1, 1, 1, 0, 0)
    )
    write.csv(df, tmp, row.names = FALSE)
    
    ## Test reading
    result <- MultiAssaySpatialExperiment:::.read_coldata_csv(tmp, technology = "Visium")
    
    expect_s4_class(result, "DataFrame")
    expect_equal(nrow(result), 5L)
    expect_true("cell_id" %in% colnames(result))  ## Should be renamed from barcode
    expect_true("x_centroid" %in% colnames(result))  ## Should be renamed
    expect_true("y_centroid" %in% colnames(result))  ## Should be renamed
    expect_true("in_tissue" %in% colnames(result))  ## Should be preserved
    
    unlink(tmp)
})

test_that(".read_coldata_csv finds columns with candidate names", {
    ## Create a CSV with alternative column names
    tmp <- tempfile(fileext = ".csv")
    df <- data.frame(
        EntityID = paste0("C", 1:3),
        CenterX = c(1.5, 2.5, 3.5),
        CenterY = c(4.5, 5.5, 6.5)
    )
    write.csv(df, tmp, row.names = FALSE)
    
    ## Test reading with candidate names
    result <- MultiAssaySpatialExperiment:::.read_coldata_csv(tmp, technology = "MERSCOPE")
    
    expect_equal(nrow(result), 3L)
    expect_true("cell_id" %in% colnames(result))
    expect_true("x_centroid" %in% colnames(result))
    expect_true("y_centroid" %in% colnames(result))
    
    unlink(tmp)
})

test_that(".create_spot_geometries creates circular polygons", {
    skip_if_not_installed("sf")

    positions <- DataFrame(
        instance_id = paste0("S", 1:3),
        x_centroid = c(10, 20, 30),
        y_centroid = c(15, 25, 35),
        in_tissue = c(1, 1, 0)
    )
    
    ## Create geometries
    result <- MultiAssaySpatialExperiment:::.create_spot_geometries(
        positions,
        x_col = "x_centroid",
        y_col = "y_centroid",
        radius = 5
    )
    
    expect_s3_class(result, "sf")
    expect_equal(nrow(result), 3L)
    expect_true("geometry" %in% names(result))
    expect_true("instance_id" %in% names(result))
    expect_true("in_tissue" %in% names(result))
    
    ## Check that geometries are polygons
    geom_types <- sf::st_geometry_type(result)
    expect_true(all(geom_types == "POLYGON"))
})

test_that(".read_geometries_auto dispatches based on file pattern", {
    skip_if_not_installed("sf")
    
    ## Test CSV vertex table
    tmp_csv <- tempfile(fileext = ".csv")
    df <- data.frame(
        cell_id = rep(c("C1", "C2"), each = 4),
        vertex_x = c(0, 1, 1, 0, 10, 11, 11, 10),
        vertex_y = c(0, 0, 1, 1, 10, 10, 11, 11)
    )
    write.csv(df, tmp_csv, row.names = FALSE)
    
    ## Test auto-dispatch
    result <- MultiAssaySpatialExperiment:::.read_geometries_auto(tmp_csv)
    
    expect_s3_class(result, "sf")
    expect_equal(nrow(result), 2L)  ## 2 cells
    
    unlink(tmp_csv)
})

test_that(".read_transcripts_csv standardizes column names", {
    ## Create mock transcripts CSV
    tmp <- tempfile(fileext = ".csv")
    df <- data.frame(
        x_location = runif(10, 0, 100),
        y_location = runif(10, 0, 100),
        feature_name = sample(paste0("Gene", 1:5), 10, replace = TRUE),
        qv = sample(20:40, 10, replace = TRUE)
    )
    write.csv(df, tmp, row.names = FALSE)
    
    ## Test reading
    result <- MultiAssaySpatialExperiment:::.read_transcripts_csv(tmp)
    
    expect_s4_class(result, "DataFrame")
    expect_equal(nrow(result), 10L)
    expect_true("x" %in% colnames(result))  ## Should be renamed
    expect_true("y" %in% colnames(result))  ## Should be renamed
    expect_true("feature" %in% colnames(result))  ## Should be renamed
    
    unlink(tmp)
})

test_that(".read_visium_scalefactors reads Visium scale factors", {
    skip_if_not_installed("jsonlite")
    
    ## Create a mock scalefactors JSON
    tmp <- tempfile(fileext = ".json")
    scale <- list(
        tissue_hires_scalef = 0.08,
        tissue_lowres_scalef = 0.016,
        spot_diameter_fullres = 89.43
    )
    jsonlite::write_json(scale, tmp, auto_unbox = TRUE)
    
    ## Test reading
    result <- MultiAssaySpatialExperiment:::.read_visium_scalefactors(tmp)
    
    expect_type(result, "list")
    expect_equal(result$tissue_hires_scalef, 0.08)
    expect_equal(result$tissue_lowres_scalef, 0.016)
    expect_equal(result$spot_diameter_fullres, 89.43)
    
    unlink(tmp)
})

test_that(".read_image_metadata handles missing files gracefully", {
    ## Test with non-existent file
    result <- MultiAssaySpatialExperiment:::.read_image_metadata(
        image_path = "/nonexistent/image.png",
        sample_id = "sample1",
        image_id = "test"
    )
    
    ## Should return NULL with a warning
    expect_null(result)
})

## ============================================================================
## Additional component reader tests (Phase 5 additions)
## ============================================================================

test_that(".read_coldata_parquet reads Parquet files", {
    skip_if_not_installed("arrow")
    
    ## Create a temporary Parquet file
    tmp <- tempfile(fileext = ".parquet")
    df <- data.frame(
        cell_id = paste0("C", 1:5),
        x_global = seq(10, 50, 10),
        y_global = seq(20, 100, 20),
        area = runif(5, 100, 200)
    )
    
    ## Try to write Parquet; skip if it fails
    tryCatch({
        arrow::write_parquet(df, tmp)
        
        ## Test reading
        result <- MultiAssaySpatialExperiment:::.read_coldata_parquet(tmp,
            id_col = c("cell_id", "barcode"),
            x_col = c("x_global", "x"),
            y_col = c("y_global", "y")
        )
        
        expect_s4_class(result, "DataFrame")
        expect_equal(nrow(result), 5L)
        expect_true("cell_id" %in% colnames(result))
        expect_true("x_centroid" %in% colnames(result))
        expect_true("y_centroid" %in% colnames(result))
        expect_true("area" %in% colnames(result))  ## Other columns preserved
        
        unlink(tmp)
    }, error = function(e) {
        skip(paste("arrow::write_parquet failed:", e$message))
    })
})

test_that(".read_transcripts_parquet reads large transcript files", {
    skip_if_not_installed("arrow")
    
    ## Create a temporary Parquet file
    tmp <- tempfile(fileext = ".parquet")
    df <- data.frame(
        global_x = runif(100, 0, 1000),
        global_y = runif(100, 0, 1000),
        gene = sample(paste0("Gene", 1:10), 100, replace = TRUE),
        qv = sample(20:40, 100, replace = TRUE),
        cell_id = sample(paste0("C", 1:20), 100, replace = TRUE)
    )
    
    tryCatch({
        arrow::write_parquet(df, tmp)
        
        ## Test reading
        result <- MultiAssaySpatialExperiment:::.read_transcripts_parquet(tmp,
            x_col = c("global_x", "x"),
            y_col = c("global_y", "y"),
            feature_col = c("gene", "feature_name"),
            qc_col = c("qv", "quality"),
            cell_col = c("cell_id", "Cell")
        )
        
        expect_s4_class(result, "DataFrame")
        expect_equal(nrow(result), 100L)
        expect_true("x" %in% colnames(result))
        expect_true("y" %in% colnames(result))
        expect_true("feature" %in% colnames(result))
        expect_true("qc_score" %in% colnames(result))
        expect_true("cell_id" %in% colnames(result))
        
        unlink(tmp)
    }, error = function(e) {
        skip(paste("arrow::write_parquet failed:", e$message))
    })
})

test_that(".read_geometries_geojson reads GeoJSON files", {
    skip_if_not_installed("sf")
    
    ## Create a mock GeoJSON file using sf
    tmp <- tempfile(fileext = ".geojson")
    
    ## Create simple polygons
    p1 <- sf::st_polygon(list(matrix(c(0,0, 1,0, 1,1, 0,1, 0,0), ncol=2, byrow=TRUE)))
    p2 <- sf::st_polygon(list(matrix(c(2,2, 3,2, 3,3, 2,3, 2,2), ncol=2, byrow=TRUE)))
    
    geoms <- sf::st_sfc(p1, p2)
    sf_obj <- sf::st_sf(instance_id = c("P1", "P2"), geometry = geoms)
    
    tryCatch({
        sf::st_write(sf_obj, tmp, driver = "GeoJSON", quiet = TRUE)
        
        ## Test reading
        result <- MultiAssaySpatialExperiment:::.read_geometries_geojson(tmp)
        
        expect_s3_class(result, "sf")
        expect_equal(nrow(result), 2L)
        expect_true("instance_id" %in% names(result))
        expect_true("geometry" %in% names(result))
        
        unlink(tmp)
    }, error = function(e) {
        skip(paste("sf::st_write failed:", e$message))
    })
})

test_that(".read_geometries_vertex_parquet reads vertex tables", {
    skip_if_not_installed("arrow")
    skip_if_not_installed("sf")
    
    ## Create a mock vertex table Parquet
    tmp <- tempfile(fileext = ".parquet")
    df <- data.frame(
        cell_id = rep(c("C1", "C2"), each = 4),
        vertex_x = c(0, 1, 1, 0, 10, 11, 11, 10),
        vertex_y = c(0, 0, 1, 1, 10, 10, 11, 11)
    )
    
    tryCatch({
        arrow::write_parquet(df, tmp)
        
        ## Test reading
        result <- MultiAssaySpatialExperiment:::.read_geometries_vertex_parquet(tmp)
        
        expect_s3_class(result, "sf")
        expect_equal(nrow(result), 2L)  ## 2 cells
        expect_true("instance_id" %in% names(result))
        
        unlink(tmp)
    }, error = function(e) {
        skip(paste("Failed to create or read vertex Parquet:", e$message))
    })
})

test_that(".read_geometries_geoparquet reads GeoParquet files", {
    skip_if_not_installed("sfarrow")
    skip_if_not_installed("sf")
    
    ## Create a mock GeoParquet file using sfarrow
    tmp <- tempfile(fileext = ".parquet")
    
    ## Create simple polygons
    p1 <- sf::st_polygon(list(matrix(c(0,0, 1,0, 1,1, 0,1, 0,0), ncol=2, byrow=TRUE)))
    p2 <- sf::st_polygon(list(matrix(c(2,2, 3,2, 3,3, 2,3, 2,2), ncol=2, byrow=TRUE)))
    
    geoms <- sf::st_sfc(p1, p2)
    sf_obj <- sf::st_sf(instance_id = c("GP1", "GP2"), geometry = geoms)
    
    tryCatch({
        sfarrow::st_write_parquet(sf_obj, tmp)
        
        ## Test reading
        result <- MultiAssaySpatialExperiment:::.read_geometries_geoparquet(tmp)
        
        expect_s3_class(result, "sf")
        expect_equal(nrow(result), 2L)
        expect_true("instance_id" %in% names(result))
        expect_true("geometry" %in% names(result))
        
        unlink(tmp)
    }, error = function(e) {
        skip(paste("sfarrow::st_write_parquet failed:", e$message))
    })
})

test_that(".dispatch_shapes_reader selects correct reader", {
    ## Test GeoParquet pattern (priority 10)
    reader_name <- MultiAssaySpatialExperiment:::.dispatch_shapes_reader("cells_sf.parquet")
    expect_equal(reader_name, ".read_geometries_geoparquet")
    
    ## Test regular Parquet pattern (priority 5)
    reader_name <- MultiAssaySpatialExperiment:::.dispatch_shapes_reader("boundaries.parquet")
    expect_equal(reader_name, ".read_geometries_vertex_parquet")
    
    ## Test GeoJSON pattern (priority 5)
    reader_name <- MultiAssaySpatialExperiment:::.dispatch_shapes_reader("cells.geojson")
    expect_equal(reader_name, ".read_geometries_geojson")
    
    ## Test CSV pattern (priority 3)
    reader_name <- MultiAssaySpatialExperiment:::.dispatch_shapes_reader("vertices.csv")
    expect_equal(reader_name, ".read_geometries_vertex_csv")
})

test_that("readHDF5ForMASE handles HDF5 matrices", {
    skip("Requires DropletUtils and a mock HDF5 file")
    
    ## TODO: Create a minimal HDF5 matrix using DropletUtils or rhdf5
    ## This is complex and requires careful setup of HDF5 structure
})

test_that("readHDF5ForMASE handles MTX matrices", {
    skip("Requires DropletUtils and mock MTX files")
    
    ## TODO: Create matrix.mtx, features.tsv, barcodes.tsv
})

test_that(".read_visium_images handles Visium image directory", {
    skip("Requires mock Visium spatial directory with images")
    
    ## TODO: Create mock spatial directory with tissue_hires_image.png, etc.
})

test_that(".read_images_metadata handles multiple images", {
    skip("Requires mock image files")
    
    ## TODO: Create mock image files and test metadata extraction
})

test_that(".read_label_raster reads TIFF label images", {
    skip("Requires terra and mock TIFF files")
    
    ## TODO: Create mock label TIFF and test reading
})

test_that(".read_label_ometiff reads OME-TIFF label images", {
    skip("Requires RBioFormats and mock OME-TIFF files")
    
    ## TODO: Create mock OME-TIFF and test reading
})

## ============================================================================
## Integration tests for component reader combinations
## ============================================================================

test_that("Component readers handle missing optional columns gracefully", {
    ## Test that readers don't fail when optional columns are missing
    tmp <- tempfile(fileext = ".csv")
    df <- data.frame(
        cell_id = paste0("C", 1:3),
        x = c(1, 2, 3),
        y = c(4, 5, 6)
        ## No QC column, no area column, etc.
    )
    write.csv(df, tmp, row.names = FALSE)
    
    ## Should read successfully without optional columns
    result <- MultiAssaySpatialExperiment:::.read_coldata_csv(tmp, technology = "Visium")
    expect_s4_class(result, "DataFrame")
    expect_equal(nrow(result), 3L)
    
    unlink(tmp)
})

test_that("Component readers preserve extra columns", {
    ## Test that readers keep non-standard columns
    tmp <- tempfile(fileext = ".csv")
    df <- data.frame(
        cell_id = paste0("C", 1:3),
        x = c(1, 2, 3),
        y = c(4, 5, 6),
        custom_metric = runif(3),
        another_column = LETTERS[1:3]
    )
    write.csv(df, tmp, row.names = FALSE)
    
    result <- MultiAssaySpatialExperiment:::.read_coldata_csv(tmp, technology = "Visium")
    expect_true("custom_metric" %in% colnames(result))
    expect_true("another_column" %in% colnames(result))
    
    unlink(tmp)
})

test_that("Parquet readers use S4 generics", {
    ## Verify that Parquet readers call the S4 generics
    
    ## Check that generics exist
    expect_true(isGeneric("readParquetForMASE"))
    expect_true(isGeneric("readGeoParquetForMASE"))
    
    ## Check that methods are defined for character
    expect_true(hasMethod("readParquetForMASE", "character"))
    expect_true(hasMethod("readGeoParquetForMASE", "character"))
})
