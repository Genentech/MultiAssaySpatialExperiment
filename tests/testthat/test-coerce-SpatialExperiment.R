test_that("as(SpatialExperiment, 'MultiAssaySpatialExperiment') works", {
    mat <- matrix(1:12, 3, 4)
    colnames(mat) <- paste0("cell", 1:4)
    spe <- SpatialExperiment::SpatialExperiment(
        list(counts = mat),
        spatialCoords = matrix(rnorm(8), 4, 2)
    )
    mase <- as(spe, "MultiAssaySpatialExperiment")
    expect_s4_class(mase, "MultiAssaySpatialExperiment")
    expect_true(validObject(mase))
    expect_equal(length(mase), 1L)
    expect_equal(names(mase), "spatial")
    expect_equal(ncol(experiments(mase)[[1L]]), 4L)
    expect_equal(length(spatialPoints(mase)), 1L)
    expect_equal(names(spatialPoints(mase)), "coordinates")
})

test_that("as(MultiAssaySpatialExperiment, 'SpatialExperiment') works with single assay", {
    mat <- matrix(1:12, 3, 4)
    colnames(mat) <- paste0("cell", 1:4)
    spe <- SpatialExperiment::SpatialExperiment(
        list(counts = mat),
        spatialCoords = matrix(rnorm(8), 4, 2)
    )
    mase <- as(spe, "MultiAssaySpatialExperiment")
    spe2 <- as(mase, "SpatialExperiment")
    expect_s4_class(spe2, "SpatialExperiment")
    expect_true(validObject(spe2))
    expect_equal(ncol(spe2), ncol(spe))
    expect_equal(nrow(spe2), nrow(spe))
})

test_that("SpatialFeatureExperiment coerces via SpatialExperiment", {
    skip_if_not_installed("SpatialFeatureExperiment")
    mat <- matrix(1:12, 3, 4)
    colnames(mat) <- paste0("cell", 1:4)
    spe <- SpatialExperiment::SpatialExperiment(
        list(counts = mat),
        spatialCoords = matrix(rnorm(8), 4, 2)
    )
    sfe <- as(spe, "SpatialFeatureExperiment")
    mase <- as(sfe, "MultiAssaySpatialExperiment")
    expect_s4_class(mase, "MultiAssaySpatialExperiment")
    expect_equal(names(mase), "spatial")
})

test_that("MASE to SpatialFeatureExperiment round-trip works", {
    skip_if_not_installed("SpatialFeatureExperiment")
    mat <- matrix(1:12, 3, 4)
    colnames(mat) <- paste0("cell", 1:4)
    spe <- SpatialExperiment::SpatialExperiment(
        list(counts = mat),
        spatialCoords = matrix(rnorm(8), 4, 2)
    )
    sfe <- as(spe, "SpatialFeatureExperiment")
    mase <- as(sfe, "MultiAssaySpatialExperiment")
    sfe2 <- as(mase, "SpatialFeatureExperiment")
    expect_s4_class(sfe2, "SpatialFeatureExperiment")
    expect_true(validObject(sfe2))
    expect_equal(ncol(sfe2), ncol(sfe))
})

test_that("empty SpatialExperiment coerces to MASE", {
    spe <- SpatialExperiment::SpatialExperiment(
        list(counts = matrix(0, 0, 0)),
        spatialCoords = matrix(numeric(0), 0, 2)
    )
    mase <- as(spe, "MultiAssaySpatialExperiment")
    expect_s4_class(mase, "MultiAssaySpatialExperiment")
    expect_true(validObject(mase))
    expect_equal(length(mase), 1L)
    expect_equal(names(mase), "spatial")
})

test_that("imgData is preserved in SPE-MASE-SPE round-trip", {
    mat <- matrix(1:12, 3, 4)
    colnames(mat) <- paste0("cell", 1:4)
    raster_img <- as.raster(matrix("#000000", 1, 1))
    idf <- DataFrame(
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
    spe2 <- as(mase, "SpatialExperiment")
    expect_false(is.null(imgData(spe2)))
    expect_equal(nrow(imgData(spe2)), nrow(imgData(spe)))
    expect_equal(imgData(spe2)$sample_id, imgData(spe)$sample_id)
})

test_that("as(MASE, 'SpatialExperiment') errors when multiple assays", {
    mat <- matrix(1:12, 3, 4, dimnames = list(paste0("G", 1:3), paste0("S", 1:4)))
    mase <- MultiAssaySpatialExperiment(
        experiments = ExperimentList(a = mat, b = mat),
        colData = DataFrame(row.names = paste0("S", 1:4)),
        sampleMap = DataFrame(
            assay = factor(c(rep("a", 4), rep("b", 4)), c("a", "b")),
            primary = rep(paste0("S", 1:4), 2),
            colname = rep(paste0("S", 1:4), 2)
        )
    )
    expect_error(as(mase, "SpatialExperiment"),
        "exactly one assay"
    )
})
