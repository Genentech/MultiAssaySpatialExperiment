### =========================================================================
### Test methods for SpatialData generics (Phase 8)
### =========================================================================

test_that("SpatialData API: images, labels, points, shapes, tables", {
    skip_if_not_installed("SpatialData")
    library(SpatialData)
    mat <- matrix(rnorm(20), nrow = 5, ncol = 4,
        dimnames = list(paste0("G", 1:5), paste0("S", 1:4)))
    img <- matrix(1:12, 3, 4)
    pts <- DataFrame(x = c(1, 2), y = c(3, 4))
    shp <- DataFrame(x = c(1, 2, 3), y = c(4, 5, 6))
    mase <- MultiAssaySpatialExperiment(
        experiments = ExperimentList(assay1 = mat),
        colData = DataFrame(row.names = paste0("S", 1:4)),
        sampleMap = DataFrame(
            assay = factor("assay1", "assay1"),
            primary = paste0("S", 1:4),
            colname = paste0("S", 1:4)
        ),
        images = RasterLayerList(brightfield = img),
        labels = RasterLayerList(cells = matrix(1:6, 2, 3)),
        points = PointsLayerList(transcripts = pts),
        shapes = ShapesLayerList(regions = shp)
    )
    expect_equal(length(images(mase)), 1L)
    expect_equal(imageNames(mase), "brightfield")
    expect_equal(image(mase, 1), img)
    expect_equal(image(mase, "brightfield"), img)
    expect_equal(length(labels(mase)), 1L)
    expect_equal(labelNames(mase), "cells")
    expect_equal(length(points(mase)), 1L)
    expect_equal(pointNames(mase), "transcripts")
    expect_equal(length(shapes(mase)), 1L)
    expect_equal(shapeNames(mase), "regions")
    expect_equal(length(tables(mase)), 1L)
    expect_equal(tableNames(mase), "assay1")
    expect_equal(SpatialData::table(mase, 1), experiments(mase)[[1L]])
})

test_that("SpatialData API: layer and element", {
    skip_if_not_installed("SpatialData")
    library(SpatialData)
    mat <- matrix(rnorm(20), nrow = 5, ncol = 4,
        dimnames = list(paste0("G", 1:5), paste0("S", 1:4)))
    img <- matrix(1:12, 3, 4)
    mase <- MultiAssaySpatialExperiment(
        experiments = ExperimentList(assay1 = mat),
        colData = DataFrame(row.names = paste0("S", 1:4)),
        sampleMap = DataFrame(
            assay = factor("assay1", "assay1"),
            primary = paste0("S", 1:4),
            colname = paste0("S", 1:4)
        ),
        images = RasterLayerList(brightfield = img)
    )
    expect_equal(layer(mase, "images"), spatialImages(mase))
    expect_equal(layer(mase, 1L), spatialImages(mase))
    expect_equal(element(mase, "images", "brightfield"), img)
    expect_equal(element(mase, "images", 1L), img)
})

test_that("SpatialData API: images<-, labels<-, etc. (batch setters)", {
    skip_if_not_installed("SpatialData")
    library(SpatialData)
    mat <- matrix(rnorm(20), nrow = 5, ncol = 4,
        dimnames = list(paste0("G", 1:5), paste0("S", 1:4)))
    mase <- MultiAssaySpatialExperiment(
        experiments = ExperimentList(assay1 = mat),
        colData = DataFrame(row.names = paste0("S", 1:4)),
        sampleMap = DataFrame(
            assay = factor("assay1", "assay1"),
            primary = paste0("S", 1:4),
            colname = paste0("S", 1:4)
        )
    )
    images(mase) <- list(newimg = matrix(1:6, 2, 3))
    expect_equal(imageNames(mase), "newimg")
    labels(mase) <- list(mask = matrix(1:8, 2, 4))
    expect_equal(labelNames(mase), "mask")
    points(mase) <- list(coords = DataFrame(x = 1:3, y = 2:4))
    expect_equal(pointNames(mase), "coords")
    shapes(mase) <- list(poly = DataFrame(a = 1:2, b = 3:4))
    expect_equal(shapeNames(mase), "poly")
})

test_that("SpatialData API: image<-, label<-, etc. (single-element setters)", {
    skip_if_not_installed("SpatialData")
    library(SpatialData)
    mat <- matrix(rnorm(20), nrow = 5, ncol = 4,
        dimnames = list(paste0("G", 1:5), paste0("S", 1:4)))
    img <- matrix(1:12, 3, 4)
    mase <- MultiAssaySpatialExperiment(
        experiments = ExperimentList(assay1 = mat),
        colData = DataFrame(row.names = paste0("S", 1:4)),
        sampleMap = DataFrame(
            assay = factor("assay1", "assay1"),
            primary = paste0("S", 1:4),
            colname = paste0("S", 1:4)
        ),
        images = RasterLayerList(a = img)
    )
    image(mase, "a") <- matrix(20:31, 3, 4)
    expect_equal(image(mase, "a")[1, 1], 20)
    image(mase, "b") <- matrix(1:6, 2, 3)
    expect_equal(imageNames(mase), c("a", "b"))
    image(mase, "a") <- NULL
    expect_equal(imageNames(mase), "b")
})
