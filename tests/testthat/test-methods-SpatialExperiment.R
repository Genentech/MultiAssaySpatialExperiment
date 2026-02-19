### =========================================================================
### Test methods for SpatialExperiment generics (Phase 6)
### =========================================================================

test_that("spatialCoords returns matrix from SPE experiment", {
    mat <- matrix(1:12, 3, 4)
    colnames(mat) <- paste0("cell", 1:4)
    scoords <- matrix(rnorm(8), 4, 2, dimnames = list(NULL, c("x", "y")))
    spe <- SpatialExperiment::SpatialExperiment(
        list(counts = mat),
        spatialCoords = scoords
    )
    mase <- as(spe, "MultiAssaySpatialExperiment")
    sc <- spatialCoords(mase)
    expect_true(is.matrix(sc))
    expect_equal(dim(sc), c(4L, 2L))
    expect_equal(colnames(sc), c("x", "y"))
    expect_equal(sc, spatialCoords(spe))
})

test_that("spatialCoords returns from points when no SPE assay", {
    pt_df <- S4Vectors::DataFrame(
        x = runif(5), y = runif(5),
        instance_id = paste0("cell", 1:5)
    )
    pts <- PointsLayerList(coordinates = pt_df)
    exp <- SummarizedExperiment::SummarizedExperiment(
        matrix(1:15, 3, 5, dimnames = list(letters[1:3], paste0("cell", 1:5)))
    )
    mase <- MultiAssaySpatialExperiment(
        experiments = ExperimentList(spatial = exp),
        colData = S4Vectors::DataFrame(row.names = "S1"),
        sampleMap = S4Vectors::DataFrame(
            assay = factor("spatial", "spatial"),
            primary = "S1",
            colname = paste0("cell", 1:5)
        ),
        points = pts
    )
    sc <- spatialCoords(mase)
    expect_true(is.matrix(sc))
    expect_equal(dim(sc), c(5L, 2L))
    expect_equal(sc[, "x"], pt_df$x)
    expect_equal(sc[, "y"], pt_df$y)
})

test_that("spatialCoordsNames returns coordinate names", {
    mat <- matrix(1:12, 3, 4)
    colnames(mat) <- paste0("cell", 1:4)
    scoords <- matrix(rnorm(8), 4, 2, dimnames = list(NULL, c("sdimx", "sdimy")))
    spe <- SpatialExperiment::SpatialExperiment(
        list(counts = mat),
        spatialCoords = scoords
    )
    mase <- as(spe, "MultiAssaySpatialExperiment")
    expect_equal(spatialCoordsNames(mase), c("sdimx", "sdimy"))
})

test_that("getImg, imgRaster, imgSource return NULL when imgData empty", {
    mat <- matrix(1:12, 3, 4)
    colnames(mat) <- paste0("cell", 1:4)
    spe <- SpatialExperiment::SpatialExperiment(
        list(counts = mat),
        spatialCoords = matrix(rnorm(8), 4, 2)
    )
    mase <- as(spe, "MultiAssaySpatialExperiment")
    expect_null(getImg(mase))
    expect_null(imgRaster(mase))
    expect_null(imgSource(mase))
})

test_that("getImg, imgRaster, imgSource work when imgData present", {
    mat <- matrix(1:12, 3, 4)
    colnames(mat) <- paste0("cell", 1:4)
    raster_img <- grDevices::as.raster(matrix("#000000", 2, 2))
    idf <- S4Vectors::DataFrame(
        sample_id = "sample01",
        image_id = "lowres",
        data = I(list(SpatialExperiment::SpatialImage(raster_img))),
        scaleFactor = 1
    )
    spe <- SpatialExperiment::SpatialExperiment(
        list(counts = mat),
        spatialCoords = matrix(rnorm(8), 4, 2),
        imgData = idf
    )
    mase <- as(spe, "MultiAssaySpatialExperiment")
    img <- getImg(mase)
    expect_s4_class(img, "VirtualSpatialImage")
    expect_true(inherits(imgRaster(mase), "raster"))
    expect_type(imgSource(mase), "character")
})

test_that("scaleFactors returns from imgData", {
    mat <- matrix(1:12, 3, 4)
    colnames(mat) <- paste0("cell", 1:4)
    raster_img <- grDevices::as.raster(matrix("#000000", 2, 2))
    idf <- S4Vectors::DataFrame(
        sample_id = "sample01",
        image_id = "lowres",
        data = I(list(SpatialExperiment::SpatialImage(raster_img))),
        scaleFactor = 2.5
    )
    spe <- SpatialExperiment::SpatialExperiment(
        list(counts = mat),
        spatialCoords = matrix(rnorm(8), 4, 2),
        imgData = idf
    )
    mase <- as(spe, "MultiAssaySpatialExperiment")
    expect_equal(scaleFactors(mase), 2.5)
})

test_that("molecules returns NULL when no molecules assay", {
    mat <- matrix(1:12, 3, 4)
    colnames(mat) <- paste0("cell", 1:4)
    spe <- SpatialExperiment::SpatialExperiment(
        list(counts = mat),
        spatialCoords = matrix(rnorm(8), 4, 2)
    )
    mase <- as(spe, "MultiAssaySpatialExperiment")
    expect_null(molecules(mase))
})

test_that("rmvImg removes imgData entry", {
    mat <- matrix(1:12, 3, 4)
    colnames(mat) <- paste0("cell", 1:4)
    raster_img <- grDevices::as.raster(matrix("#000000", 2, 2))
    idf <- S4Vectors::DataFrame(
        sample_id = "sample01",
        image_id = "lowres",
        data = I(list(SpatialExperiment::SpatialImage(raster_img))),
        scaleFactor = 1
    )
    spe <- SpatialExperiment::SpatialExperiment(
        list(counts = mat),
        spatialCoords = matrix(rnorm(8), 4, 2),
        imgData = idf
    )
    mase <- as(spe, "MultiAssaySpatialExperiment")
    expect_equal(nrow(imgData(mase)), 1L)
    mase2 <- rmvImg(mase, sample_id = "sample01", image_id = "lowres")
    expect_equal(nrow(imgData(mase2)), 0L)
})

test_that("rotateImg and mirrorImg modify imgData in place", {
    mat <- matrix(1:12, 3, 4)
    colnames(mat) <- paste0("cell", 1:4)
    raster_img <- grDevices::as.raster(matrix(
        c("#FF0000", "#00FF00", "#0000FF", "#FFFF00"), 2, 2))
    idf <- S4Vectors::DataFrame(
        sample_id = "sample01",
        image_id = "lowres",
        data = I(list(SpatialExperiment::SpatialImage(raster_img))),
        scaleFactor = 1
    )
    spe <- SpatialExperiment::SpatialExperiment(
        list(counts = mat),
        spatialCoords = matrix(rnorm(8), 4, 2),
        imgData = idf
    )
    mase <- as(spe, "MultiAssaySpatialExperiment")
    expect_s4_class(rotateImg(mase, degrees = 90), "MultiAssaySpatialExperiment")
    expect_s4_class(mirrorImg(mase, axis = "h"), "MultiAssaySpatialExperiment")
    expect_false(is.null(imgData(mase)))
})
