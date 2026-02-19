test_that("MultiAssaySpatialExperiment empty constructor", {
    mase <- MultiAssaySpatialExperiment()
    expect_s4_class(mase, "MultiAssaySpatialExperiment")
    expect_s4_class(mase, "MultiAssayExperiment")
    expect_true(validObject(mase))
    expect_length(experiments(mase), 0L)
    expect_length(spatialImages(mase), 0L)
    expect_length(spatialPoints(mase), 0L)
    expect_null(imgData(mase))
    expect_null(spatialMap(mase))
})

test_that("MultiAssaySpatialExperiment with MAE-style input", {
    mat <- matrix(rnorm(20), nrow = 5, ncol = 4,
        dimnames = list(paste0("G", 1:5), paste0("S", 1:4)))
    expList <- ExperimentList(assay1 = mat)
    cd <- DataFrame(row.names = paste0("S", 1:4))
    sm <- DataFrame(
        assay = factor("assay1", "assay1"),
        primary = paste0("S", 1:4),
        colname = paste0("S", 1:4)
    )
    mase <- MultiAssaySpatialExperiment(
        experiments = expList,
        colData = cd,
        sampleMap = sm
    )
    expect_true(validObject(mase))
    expect_equal(length(experiments(mase)), 1L)
    expect_equal(names(mase), "assay1")
})

test_that("MultiAssaySpatialExperiment with spatial elements", {
    mat <- matrix(rnorm(20), nrow = 5, ncol = 4,
        dimnames = list(paste0("G", 1:5), paste0("S", 1:4)))
    img <- matrix(1:12, 3, 4)
    pts <- DataFrame(x = c(1, 2), y = c(3, 4))
    expList <- ExperimentList(assay1 = mat)
    cd <- DataFrame(row.names = paste0("S", 1:4))
    sm <- DataFrame(
        assay = factor("assay1", "assay1"),
        primary = paste0("S", 1:4),
        colname = paste0("S", 1:4)
    )
    mase <- MultiAssaySpatialExperiment(
        experiments = expList,
        colData = cd,
        sampleMap = sm,
        images = RasterLayerList(brightfield = img),
        points = PointsLayerList(transcripts = pts)
    )
    expect_true(validObject(mase))
    expect_equal(names(spatialImages(mase)), "brightfield")
    expect_equal(names(spatialPoints(mase)), "transcripts")
    expect_equal(spatialImages(mase)[[1L]], img)
    expect_equal(spatialPoints(mase)[[1L]], pts)
})

test_that("show method runs without error", {
    mase <- MultiAssaySpatialExperiment()
    expect_output(show(mase), "MultiAssaySpatialExperiment")
    expect_output(show(mase), "Spatial elements")
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Phase 3: Slot setters
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

test_that("spatialImages<- replaces slot correctly", {
    mat <- matrix(rnorm(20), nrow = 5, ncol = 4,
        dimnames = list(paste0("G", 1:5), paste0("S", 1:4)))
    img <- matrix(1:12, 3, 4)
    expList <- ExperimentList(assay1 = mat)
    cd <- DataFrame(row.names = paste0("S", 1:4))
    sm <- DataFrame(
        assay = factor("assay1", "assay1"),
        primary = paste0("S", 1:4),
        colname = paste0("S", 1:4)
    )
    mase <- MultiAssaySpatialExperiment(
        experiments = expList,
        colData = cd,
        sampleMap = sm,
        images = RasterLayerList(brightfield = img)
    )
    new_img <- matrix(20:31, 3, 4)
    spatialImages(mase) <- RasterLayerList(fluorescent = new_img)
    expect_equal(names(spatialImages(mase)), "fluorescent")
    expect_equal(spatialImages(mase)[[1L]], new_img)
    expect_true(validObject(mase))
})

test_that("spatialLabels<-, spatialPoints<-, spatialShapes<- work", {
    mat <- matrix(rnorm(20), nrow = 5, ncol = 4,
        dimnames = list(paste0("G", 1:5), paste0("S", 1:4)))
    pts <- DataFrame(x = c(1, 2), y = c(3, 4))
    shp <- DataFrame(x = c(1, 2, 3), y = c(4, 5, 6))
    expList <- ExperimentList(assay1 = mat)
    cd <- DataFrame(row.names = paste0("S", 1:4))
    sm <- DataFrame(
        assay = factor("assay1", "assay1"),
        primary = paste0("S", 1:4),
        colname = paste0("S", 1:4)
    )
    mase <- MultiAssaySpatialExperiment(
        experiments = expList,
        colData = cd,
        sampleMap = sm,
        points = PointsLayerList(transcripts = pts)
    )
    lbl <- matrix(1:6, 2, 3)
    spatialLabels(mase) <- RasterLayerList(cells = lbl)
    expect_equal(names(spatialLabels(mase)), "cells")
    expect_equal(spatialLabels(mase)[[1L]], lbl)

    spatialPoints(mase) <- PointsLayerList(spots = pts)
    expect_equal(names(spatialPoints(mase)), "spots")

    spatialShapes(mase) <- ShapesLayerList(regions = shp)
    expect_equal(names(spatialShapes(mase)), "regions")
    expect_true(validObject(mase))
})

test_that("imgData<- and spatialMap<- work", {
    mase <- MultiAssaySpatialExperiment()
    idf <- DataFrame(sample_id = "S1", data_id = "img1")
    imgData(mase) <- idf
    expect_equal(imgData(mase), idf)
    expect_true(validObject(mase))

    smap <- DataFrame(
        assay = factor("a", "a"),
        colname = "c",
        region = "r1",
        instance_id = 1L
    )
    spatialMap(mase) <- smap
    expect_equal(spatialMap(mase), smap)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Phase 3: Single-element accessors
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

test_that("single-element accessors return NULL for non-matching character index", {
    mase <- MultiAssaySpatialExperiment()
    expect_null(spatialImage(mase, "nonexistent"))
})

test_that("spatialImage, spatialLabel, spatialPoint, spatialShape work", {
    mat <- matrix(rnorm(20), nrow = 5, ncol = 4,
        dimnames = list(paste0("G", 1:5), paste0("S", 1:4)))
    img <- matrix(1:12, 3, 4)
    lbl <- matrix(1:6, 2, 3)
    pts <- DataFrame(x = c(1, 2), y = c(3, 4))
    shp <- DataFrame(x = c(1, 2, 3), y = c(4, 5, 6))
    expList <- ExperimentList(assay1 = mat)
    cd <- DataFrame(row.names = paste0("S", 1:4))
    sm <- DataFrame(
        assay = factor("assay1", "assay1"),
        primary = paste0("S", 1:4),
        colname = paste0("S", 1:4)
    )
    mase <- MultiAssaySpatialExperiment(
        experiments = expList,
        colData = cd,
        sampleMap = sm,
        images = RasterLayerList(brightfield = img),
        labels = RasterLayerList(cells = lbl),
        points = PointsLayerList(transcripts = pts),
        shapes = ShapesLayerList(regions = shp)
    )
    expect_equal(spatialImage(mase, 1L), img)
    expect_equal(spatialImage(mase, "brightfield"), img)
    expect_equal(spatialLabel(mase, "cells"), lbl)
    expect_equal(spatialPoint(mase, "transcripts"), pts)
    expect_equal(spatialShape(mase, "regions"), shp)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Phase 3: Single-element setters
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

test_that("spatialImage<-, spatialLabel<-, spatialPoint<-, spatialShape<- work", {
    mat <- matrix(rnorm(20), nrow = 5, ncol = 4,
        dimnames = list(paste0("G", 1:5), paste0("S", 1:4)))
    img <- matrix(1:12, 3, 4)
    pts <- DataFrame(x = c(1, 2), y = c(3, 4))
    expList <- ExperimentList(assay1 = mat)
    cd <- DataFrame(row.names = paste0("S", 1:4))
    sm <- DataFrame(
        assay = factor("assay1", "assay1"),
        primary = paste0("S", 1:4),
        colname = paste0("S", 1:4)
    )
    mase <- MultiAssaySpatialExperiment(
        experiments = expList,
        colData = cd,
        sampleMap = sm,
        images = RasterLayerList(brightfield = img),
        points = PointsLayerList(transcripts = pts)
    )
    new_img <- matrix(99:110, 3, 4)
    spatialImage(mase, "brightfield") <- new_img
    expect_equal(spatialImage(mase, "brightfield"), new_img)

    new_pts <- DataFrame(x = c(10, 20), y = c(30, 40))
    spatialPoint(mase, "transcripts") <- new_pts
    expect_equal(spatialPoint(mase, "transcripts"), new_pts)

    lbl <- matrix(1:6, 2, 3)
    spatialLabel(mase, "cells") <- lbl
    expect_equal(spatialLabel(mase, "cells"), lbl)

    shp <- DataFrame(x = c(1, 2), y = c(3, 4))
    spatialShape(mase, "regions") <- shp
    expect_equal(spatialShape(mase, "regions"), shp)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Phase 3: Name accessors
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

test_that("spatialImageNames, spatialLabelNames, spatialPointNames, spatialShapeNames work", {
    mat <- matrix(rnorm(20), nrow = 5, ncol = 4,
        dimnames = list(paste0("G", 1:5), paste0("S", 1:4)))
    img <- matrix(1:12, 3, 4)
    lbl <- matrix(1:6, 2, 3)
    pts <- DataFrame(x = c(1, 2), y = c(3, 4))
    shp <- DataFrame(x = c(1, 2, 3), y = c(4, 5, 6))
    expList <- ExperimentList(assay1 = mat)
    cd <- DataFrame(row.names = paste0("S", 1:4))
    sm <- DataFrame(
        assay = factor("assay1", "assay1"),
        primary = paste0("S", 1:4),
        colname = paste0("S", 1:4)
    )
    mase <- MultiAssaySpatialExperiment(
        experiments = expList,
        colData = cd,
        sampleMap = sm,
        images = RasterLayerList(brightfield = img),
        labels = RasterLayerList(cells = lbl),
        points = PointsLayerList(transcripts = pts),
        shapes = ShapesLayerList(regions = shp)
    )
    expect_equal(spatialImageNames(mase), "brightfield")
    expect_equal(spatialLabelNames(mase), "cells")
    expect_equal(spatialPointNames(mase), "transcripts")
    expect_equal(spatialShapeNames(mase), "regions")
})
