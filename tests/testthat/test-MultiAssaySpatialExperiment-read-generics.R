test_that("S4 generics are defined", {
    ## Test that all 5 generics exist
    expect_true(isGeneric("readParquetForMASE"))
    expect_true(isGeneric("readGeoParquetForMASE"))
    expect_true(isGeneric("readHDF5ForMASE"))
    expect_true(isGeneric("readGeoJSONForMASE"))
    expect_true(isGeneric("readCSVForMASE"))
})

test_that("readCSVForMASE works with mock data", {
    ## Create a temporary CSV file
    tmp <- tempfile(fileext = ".csv")
    df <- data.frame(
        cell_id = paste0("C", 1:5),
        x_centroid = runif(5, 0, 100),
        y_centroid = runif(5, 0, 100),
        gene_count = rpois(5, 100)
    )
    write.csv(df, tmp, row.names = FALSE)
    
    ## Test reading
    result <- readCSVForMASE(tmp)
    
    expect_s4_class(result, "DataFrame")
    expect_equal(nrow(result), 5L)
    expect_equal(ncol(result), 4L)
    expect_true(all(c("cell_id", "x_centroid", "y_centroid", "gene_count") %in% colnames(result)))
    
    unlink(tmp)
})

test_that("readCSVForMASE handles missing file gracefully", {
    expect_error(readCSVForMASE("/nonexistent/path/file.csv"))
})

test_that("readParquetForMASE requires arrow package", {
    skip_if_not_installed("arrow")
    
    ## Skip if arrow can't write parquet (some systems have issues)
    skip_on_cran()
    
    ## Create a temporary parquet file with explicit write
    tmp <- tempfile(fileext = ".parquet")
    df <- data.frame(
        cell_id = paste0("C", 1:10),
        x = runif(10, 0, 100),
        y = runif(10, 0, 100),
        stringsAsFactors = FALSE
    )
    
    ## Try to write and read - skip if fails (environment issue)
    skip_if(
        !tryCatch({
            arrow::write_parquet(df, tmp)
            result <- arrow::read_parquet(tmp)
            nrow(result) == 10 && ncol(result) == 3
        }, error = function(e) FALSE),
        "Arrow parquet write/read not working in this environment"
    )
    
    ## Test reading with our generic
    result <- readParquetForMASE(tmp)
    
    expect_s4_class(result, "DataFrame")
    expect_equal(nrow(result), 10L)
    expect_true(ncol(result) >= 3L)
    expect_true(all(c("cell_id", "x", "y") %in% colnames(result)))
    
    unlink(tmp)
})

test_that("readGeoJSONForMASE requires sf package", {
    skip_if_not_installed("sf")
    
    ## Create a simple GeoJSON using sf directly (more reliable)
    tmp <- tempfile(fileext = ".geojson")
    pt <- sf::st_point(c(10, 20))
    sfc <- sf::st_sfc(pt)
    sf_obj <- sf::st_sf(id = "C1", geometry = sfc)
    sf::st_write(sf_obj, tmp, quiet = TRUE, delete_dsn = TRUE)
    
    ## Test reading
    result <- readGeoJSONForMASE(tmp, quiet = TRUE)
    
    expect_s3_class(result, "sf")
    expect_equal(nrow(result), 1L)
    expect_true("geometry" %in% names(result))
    
    unlink(tmp)
})

test_that("readHDF5ForMASE requires DropletUtils for 10x format", {
    skip_if_not_installed("DropletUtils")
    skip("Requires real HDF5 test file")
})
