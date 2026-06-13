## Register SpatialData coercion when the optional package loads.
.has_spatialdata_coercion <- function() {
    if (!requireNamespace("SpatialData", quietly = TRUE))
        return(FALSE)
    if (!isNamespaceLoaded("SpatialData"))
        loadNamespace("SpatialData")
    tryCatch({
        methods::getMethod("coerce", c("SpatialData", "MultiAssaySpatialExperiment"))
        TRUE
    }, error = function(e) FALSE)
}

test_that("as(SpatialData, 'MultiAssaySpatialExperiment') works when SpatialData available", {
    skip_if_not_installed("SpatialData")
    skip_if(!.has_spatialdata_coercion(),
            "SpatialData coercion not registered")
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
    skip_if(!.has_spatialdata_coercion(),
            "SpatialData coercion not registered")
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
    skip_if(!.has_spatialdata_coercion(),
            "SpatialData coercion not registered")
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
    skip_if(!.has_spatialdata_coercion(),
            "SpatialData coercion not registered")
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
    skip_if(!.has_spatialdata_coercion(),
            "SpatialData coercion not registered")
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

test_that("SpatialData→MASE coercion includes element_type in spatialMap", {
    skip_if_not_installed("SpatialData")
    skip_if(!.has_spatialdata_coercion(),
            "SpatialData coercion not registered")
    library(SingleCellExperiment)
    
    ## Create SCE with spatialdata_attrs metadata
    sce <- SingleCellExperiment(
        list(counts = matrix(1:12, 3, 4)),
        colData = DataFrame(
            rk = rep("cells", 4),
            instance_id = paste0("cell", 1:4),
            row.names = paste0("cell", 1:4)
        )
    )
    int_metadata(sce)[["spatialdata_attrs"]] <- list(
        region = "cells",
        region_key = "rk",
        instance_key = "instance_id"
    )
    
    ## Create SpatialData with points$cells
    sd <- SpatialData::SpatialData(
        images = list(),
        labels = list(),
        points = list(
            cells = SpatialData::PointFrame(data = data.frame(
                x = 1:4, y = 2:5, instance_id = paste0("cell", 1:4)
            ))
        ),
        shapes = list(),
        tables = list(assay1 = sce)
    )
    
    ## Convert to MASE
    mase <- as(sd, "MultiAssaySpatialExperiment")
    expect_s4_class(mase, "MultiAssaySpatialExperiment")
    
    ## Check spatialMap has element_type column
    smap <- spatialMap(mase)
    expect_true(!is.null(smap))
    expect_true("element_type" %in% colnames(smap))
    
    ## Check element_type is "points" (region found in points slot)
    expect_equal(unique(smap$element_type), "points")
    expect_equal(unique(smap$region), "cells")
})

test_that("SpatialData→MASE coercion sets element_type=shapes for shape regions", {
    skip_if_not_installed("SpatialData")
    skip_if(!.has_spatialdata_coercion(),
            "SpatialData coercion not registered")
    library(SingleCellExperiment)
    
    ## Create SCE with spatialdata_attrs metadata pointing to shapes
    sce <- SingleCellExperiment(
        list(counts = matrix(1:6, 2, 3)),
        colData = DataFrame(
            rk = rep("nuclei", 3),
            instance_id = paste0("N", 1:3),
            row.names = paste0("N", 1:3)
        )
    )
    int_metadata(sce)[["spatialdata_attrs"]] <- list(
        region = "nuclei",
        region_key = "rk",
        instance_key = "instance_id"
    )
    
    ## Create SpatialData with shapes$nuclei (not points)
    sd <- SpatialData::SpatialData(
        images = list(),
        labels = list(),
        points = list(),
        shapes = list(
            nuclei = SpatialData::ShapeFrame(data = data.frame(
                id = paste0("N", 1:3),
                area = c(100, 150, 200)
            ))
        ),
        tables = list(assay1 = sce)
    )
    
    ## Convert to MASE
    mase <- as(sd, "MultiAssaySpatialExperiment")
    
    ## Check spatialMap has element_type = "shapes"
    smap <- spatialMap(mase)
    expect_true(!is.null(smap))
    expect_true("element_type" %in% colnames(smap))
    expect_equal(unique(smap$element_type), "shapes")
    expect_equal(unique(smap$region), "nuclei")
})

test_that("MASE→SD→MASE round-trip preserves element_type", {
    skip_if_not_installed("SpatialData")
    skip_if(!.has_spatialdata_coercion(),
            "SpatialData coercion not registered")
    library(SingleCellExperiment)
    
    ## Create MASE with spatialMap including element_type
    sce <- SingleCellExperiment(
        list(counts = matrix(1:12, 3, 4)),
        colData = DataFrame(row.names = paste0("cell", 1:4))
    )
    mase1 <- MultiAssaySpatialExperiment(
        experiments = ExperimentList(assay1 = sce),
        colData = DataFrame(row.names = "sample01"),
        sampleMap = DataFrame(
            assay = factor("assay1", "assay1"),
            primary = "sample01",
            colname = paste0("cell", 1:4)
        ),
        points = PointsLayerList(cells = DataFrame(
            x = c(1, 2, 3, 4),
            y = c(5, 6, 7, 8),
            instance_id = paste0("cell", 1:4)
        )),
        spatialMap = DataFrame(
            assay = factor(rep("assay1", 4), levels = "assay1"),
            colname = paste0("cell", 1:4),
            element_type = rep("points", 4),
            region = rep("cells", 4),
            instance_id = paste0("cell", 1:4)
        )
    )
    
    ## Convert to SpatialData and back
    sd <- as(mase1, "SpatialData")
    mase2 <- as(sd, "MultiAssaySpatialExperiment")
    
    ## Check element_type preserved
    smap1 <- spatialMap(mase1)
    smap2 <- spatialMap(mase2)
    
    expect_true(!is.null(smap2))
    expect_true("element_type" %in% colnames(smap2))
    expect_equal(smap2$element_type, smap1$element_type)
    expect_equal(smap2$region, smap1$region)
})
