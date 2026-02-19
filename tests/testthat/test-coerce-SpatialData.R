test_that("as(SpatialData, 'MultiAssaySpatialExperiment') works when SpatialData available", {
    skip_if_not_installed("SpatialData")
    library(SingleCellExperiment)
    sce <- SingleCellExperiment(
        list(counts = matrix(1:12, 3, 4)),
        colData = DataFrame(row.names = paste0("cell", 1:4))
    )
    sd <- SpatialData::SpatialData(
        images = list(),
        labels = list(),
        points = list(
            pts1 = SpatialData::PointFrame(data = data.frame(
                x = c(1, 2, 3, 4),
                y = c(1, 2, 3, 4),
                instance_id = 1:4
            ))
        ),
        shapes = list(
            shp1 = SpatialData::ShapeFrame(data = data.frame(
                id = 1:2,
                value = c("a", "b")
            ))
        ),
        tables = list(assay1 = sce)
    )
    mase <- as(sd, "MultiAssaySpatialExperiment")
    expect_s4_class(mase, "MultiAssaySpatialExperiment")
    expect_true(validObject(mase))
    expect_equal(length(mase), 1L)
    expect_equal(names(mase), "assay1")
    expect_equal(ncol(experiments(mase)[[1L]]), 4L)
    expect_equal(spatialPointNames(mase), "pts1")
    expect_equal(spatialShapeNames(mase), "shp1")
})

test_that("as(MultiAssaySpatialExperiment, 'SpatialData') works", {
    skip_if_not_installed("SpatialData")
    library(SingleCellExperiment)
    sce <- SingleCellExperiment(
        list(counts = matrix(1:12, 3, 4)),
        colData = DataFrame(row.names = paste0("cell", 1:4))
    )
    mase <- MultiAssaySpatialExperiment(
        experiments = ExperimentList(assay1 = sce),
        colData = DataFrame(row.names = "sample01"),
        sampleMap = DataFrame(
            assay = factor("assay1", "assay1"),
            primary = "sample01",
            colname = paste0("cell", 1:4)
        ),
        points = PointsLayerList(coords = DataFrame(
            x = c(1, 2, 3, 4),
            y = c(1, 2, 3, 4),
            instance_id = paste0("cell", 1:4)
        ))
    )
    sd <- as(mase, "SpatialData")
    expect_s4_class(sd, "SpatialData")
    expect_equal(length(SpatialData::tables(sd)), 1L)
    expect_equal(SpatialData::tableNames(sd), "assay1")
    expect_equal(SpatialData::pointNames(sd), "coords")
})

test_that("SD-MASE-SD round-trip preserves tables and points", {
    skip_if_not_installed("SpatialData")
    library(SingleCellExperiment)
    sce <- SingleCellExperiment(
        list(counts = matrix(1:12, 3, 4)),
        colData = DataFrame(row.names = paste0("cell", 1:4))
    )
    sd1 <- SpatialData::SpatialData(
        images = list(),
        labels = list(),
        points = list(
            p1 = SpatialData::PointFrame(data = data.frame(
                x = 1:4, y = 2:5, id = 1:4
            ))
        ),
        shapes = list(),
        tables = list(t1 = sce)
    )
    mase <- as(sd1, "MultiAssaySpatialExperiment")
    sd2 <- as(mase, "SpatialData")
    expect_equal(SpatialData::tableNames(sd2), "t1")
    expect_equal(ncol(SpatialData::table(sd2, "t1")), 4L)
})

test_that("multi-table SpatialData coerces to multi-assay MASE", {
    skip_if_not_installed("SpatialData")
    library(SingleCellExperiment)
    sce1 <- SingleCellExperiment(
        list(counts = matrix(1:6, 2, 3)),
        colData = DataFrame(row.names = paste0("A", 1:3))
    )
    sce2 <- SingleCellExperiment(
        list(counts = matrix(7:12, 2, 3)),
        colData = DataFrame(row.names = paste0("B", 1:3))
    )
    sd <- SpatialData::SpatialData(
        images = list(), labels = list(), points = list(), shapes = list(),
        tables = list(assay1 = sce1, assay2 = sce2)
    )
    mase <- as(sd, "MultiAssaySpatialExperiment")
    expect_s4_class(mase, "MultiAssaySpatialExperiment")
    expect_equal(length(mase), 2L)
    expect_equal(names(mase), c("assay1", "assay2"))
    expect_equal(ncol(experiments(mase)[[1L]]), 3L)
    expect_equal(ncol(experiments(mase)[[2L]]), 3L)
    expect_equal(nrow(sampleMap(mase)), 6L)
})

test_that("as(SpatialData, 'MultiAssaySpatialExperiment') errors when no tables", {
    skip_if_not_installed("SpatialData")
    sd <- SpatialData::SpatialData(
        images = list(),
        labels = list(),
        points = list(),
        shapes = list(),
        tables = list()
    )
    expect_error(as(sd, "MultiAssaySpatialExperiment"),
        "at least one table"
    )
})
