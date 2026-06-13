test_that("[ subsets MultiAssaySpatialExperiment with spatial propagation", {
    mat1 <- matrix(rnorm(20), 5, 4,
        dimnames = list(paste0("G", 1:5), paste0("S", 1:4)))
    mat2 <- matrix(rnorm(15), 3, 5,
        dimnames = list(paste0("G", 1:3), paste0("S", c(1, 2, 3, 4, 5))))
    expList <- ExperimentList(assay1 = mat1, assay2 = mat2)
    cd <- DataFrame(row.names = paste0("S", 1:5))
    sm <- DataFrame(
        assay = factor(c(rep("assay1", 4), rep("assay2", 5)),
            levels = c("assay1", "assay2")),
        primary = c(paste0("S", 1:4), paste0("S", 1:5)),
        colname = c(paste0("S", 1:4), paste0("S", 1:5))
    )
    pts <- DataFrame(
        x = c(1, 2, 3), y = c(1, 2, 3),
        instance_id = paste0("pt", 1:3)
    )
    spmap <- DataFrame(
        assay = c("assay1", "assay1", "assay2"),
        colname = c("S1", "S2", "S3"),
        element_type = "points",
        region = "coords",
        instance_id = paste0("pt", 1:3)
    )
    imgdf <- DataFrame(sample_id = paste0("S", 1:3))
    mase <- MultiAssaySpatialExperiment(
        experiments = expList,
        colData = cd,
        sampleMap = sm,
        points = PointsLayerList(coords = pts),
        images = RasterLayerList(assay1 = matrix(1:4, 2, 2)),
        spatialMap = spmap,
        imgData = imgdf
    )

    ## subset by j (colData)
    out <- mase[, c("S1", "S2")]
    expect_s4_class(out, "MultiAssaySpatialExperiment")
    expect_equal(rownames(colData(out)), c("S1", "S2"))
    expect_equal(nrow(spatialMap(out)), 2L)
    expect_equal(nrow(imgData(out)), 2L)
    expect_equal(nrow(spatialPoints(out)[["coords"]]), 2L)
    expect_equal(spatialPoints(out)[["coords"]][["instance_id"]], c("pt1", "pt2"))

    ## subset by k (assay)
    out2 <- mase[, , "assay1"]
    expect_equal(names(out2), "assay1")
    expect_equal(names(spatialImages(out2)), "assay1")
    expect_equal(unique(spatialMap(out2)[["assay"]]), "assay1")
})

test_that("subsetByAssay propagates to spatial element lists", {
    mase <- MultiAssaySpatialExperiment(
        experiments = ExperimentList(
            a = matrix(1:6, 2, 3, dimnames = list(c("R1","R2"), c("C1","C2","C3"))),
            b = matrix(1:4, 2, 2, dimnames = list(c("R1","R2"), c("C1","C2")))
        ),
        colData = DataFrame(row.names = c("C1", "C2", "C3")),
        sampleMap = DataFrame(
            assay = factor(c(rep("a", 3), rep("b", 2))),
            primary = c("C1", "C2", "C3", "C1", "C2"),
            colname = c("C1", "C2", "C3", "C1", "C2")
        ),
        images = RasterLayerList(a = matrix(1, 1, 1), b = matrix(2, 1, 1))
    )
    out <- subsetByAssay(mase, "a")
    expect_equal(names(out), "a")
    expect_equal(names(spatialImages(out)), "a")
})

test_that("c combines two MultiAssaySpatialExperiment", {
    m1 <- MultiAssaySpatialExperiment(
        experiments = ExperimentList(a = matrix(1:4, 2, 2,
            dimnames = list(c("R1","R2"), c("S1","S2")))),
        colData = DataFrame(row.names = c("S1", "S2")),
        sampleMap = DataFrame(
            assay = factor("a", "a"),
            primary = c("S1", "S2"),
            colname = c("S1", "S2")
        ),
        images = RasterLayerList(a = matrix(1, 1, 1))
    )
    m2 <- MultiAssaySpatialExperiment(
        experiments = ExperimentList(b = matrix(5:8, 2, 2,
            dimnames = list(c("R1","R2"), c("S3","S4")))),
        colData = DataFrame(row.names = c("S3", "S4")),
        sampleMap = DataFrame(
            assay = factor("b", "b"),
            primary = c("S3", "S4"),
            colname = c("S3", "S4")
        ),
        images = RasterLayerList(b = matrix(2, 1, 1))
    )
    combined <- c(m1, m2)
    expect_s4_class(combined, "MultiAssaySpatialExperiment")
    expect_equal(names(combined), c("a", "b"))
    expect_equal(nrow(colData(combined)), 4L)
    expect_equal(names(spatialImages(combined)), c("a", "b"))
})

test_that("subsetByBoundingBox filters points and propagates to assays", {
    pts <- DataFrame(x = c(1, 2, 3, 4, 5), y = c(1, 2, 3, 4, 5),
        instance_id = paste0("S", 1:5))
    mase <- MultiAssaySpatialExperiment(
        experiments = ExperimentList(
            assay1 = matrix(rnorm(25), 5, 5, dimnames = list(paste0("G", 1:5), paste0("S", 1:5)))
        ),
        colData = DataFrame(row.names = paste0("S", 1:5)),
        sampleMap = DataFrame(
            assay = factor("assay1"),
            primary = paste0("S", 1:5),
            colname = paste0("S", 1:5)
        ),
        points = PointsLayerList(coords = pts),
        spatialMap = DataFrame(
            assay = factor("assay1"),
            colname = paste0("S", 1:5),
            element_type = "points",
            region = "coords",
            instance_id = paste0("S", 1:5)
        )
    )
    out <- subsetByBoundingBox(mase, xmin = 1.5, xmax = 4.5, ymin = 1.5, ymax = 4.5)
    expect_s4_class(out, "MultiAssaySpatialExperiment")
    expect_equal(nrow(spatialPoints(out)[["coords"]]), 3L)
    expect_equal(ncol(experiments(out)[["assay1"]]), 3L)
    expect_equal(colnames(experiments(out)[["assay1"]]), c("S2", "S3", "S4"))
})

test_that("subsetByBoundingBox crops matrix rasters in spatialImages and spatialLabels", {
    pts <- DataFrame(x = c(1, 2, 3, 4, 5), y = c(1, 2, 3, 4, 5),
        instance_id = paste0("S", 1:5))
    img <- matrix(seq_len(16), 4, 4)
    mase <- MultiAssaySpatialExperiment(
        experiments = ExperimentList(
            assay1 = matrix(rnorm(25), 5, 5,
                dimnames = list(paste0("G", 1:5), paste0("S", 1:5)))
        ),
        colData = DataFrame(row.names = paste0("S", 1:5)),
        sampleMap = DataFrame(
            assay = factor("assay1"),
            primary = paste0("S", 1:5),
            colname = paste0("S", 1:5)
        ),
        points = PointsLayerList(coords = pts),
        images = RasterLayerList(tissue = img),
        labels = RasterLayerList(segmentation = img),
        spatialMap = DataFrame(
            assay = factor("assay1"),
            colname = paste0("S", 1:5),
            element_type = "points",
            region = "coords",
            instance_id = paste0("S", 1:5)
        )
    )
    out <- subsetByBoundingBox(mase, xmin = 2, xmax = 3, ymin = 2, ymax = 3)
    expect_equal(dim(spatialImages(out)[["tissue"]]), c(2L, 2L))
    expect_equal(spatialImages(out)[["tissue"]], img[2:3, 2:3, drop = FALSE])
    expect_equal(dim(spatialLabels(out)[["segmentation"]]), c(2L, 2L))
})

test_that("subsetByPolygon crops rasters to polygon bounding box", {
    skip_if_not_installed("sf")
    pts <- DataFrame(x = c(1, 2, 3), y = c(1, 2, 3), instance_id = c("S1", "S2", "S3"))
    img <- matrix(seq_len(9), 3, 3)
    mase <- MultiAssaySpatialExperiment(
        experiments = ExperimentList(
            assay1 = matrix(1:9, 3, 3, dimnames = list(paste0("G", 1:3), paste0("S", 1:3)))
        ),
        colData = DataFrame(row.names = paste0("S", 1:3)),
        sampleMap = DataFrame(
            assay = factor("assay1"),
            primary = paste0("S", 1:3),
            colname = paste0("S", 1:3)
        ),
        points = PointsLayerList(coords = pts),
        images = RasterLayerList(tissue = img),
        spatialMap = DataFrame(
            assay = factor("assay1"),
            colname = paste0("S", 1:3),
            element_type = "points",
            region = "coords",
            instance_id = paste0("S", 1:3)
        )
    )
    poly <- sf::st_polygon(list(matrix(c(2, 2, 3, 2, 3, 3, 2, 3, 2, 2),
                                     ncol = 2, byrow = TRUE)))
    out <- subsetByPolygon(mase, poly)
    expect_equal(dim(spatialImages(out)[["tissue"]]), c(2L, 2L))
    expect_equal(spatialImages(out)[["tissue"]], img[2:3, 2:3, drop = FALSE])
})

test_that("subsetByColumn propagates to spatial element lists", {
    pts <- DataFrame(
        x = 1:6, y = 1:6,
        instance_id = paste0("spot", 1:6)
    )
    ex <- SummarizedExperiment(
        assays = list(counts = matrix(1:30, 5, 6,
            dimnames = list(paste0("G", 1:5), paste0("spot", 1:6)))),
        colData = DataFrame(
            tissue = c("tumor", "tumor", "tumor", "normal", "normal", "normal"),
            nCount = c(100, 200, 150, 80, 120, 90)
        )
    )
    mase <- MultiAssaySpatialExperiment(
        experiments = ExperimentList(rna = ex),
        colData = DataFrame(row.names = c("S1", "S2")),
        sampleMap = DataFrame(
            assay = factor(rep("rna", 6), "rna"),
            primary = c(rep("S1", 3), rep("S2", 3)),
            colname = paste0("spot", 1:6)
        ),
        points = PointsLayerList(tissue1 = pts),
        spatialMap = DataFrame(
            assay = factor(rep("rna", 6), "rna"),
            colname = paste0("spot", 1:6),
            element_type = "points",
            region = "tissue1",
            instance_id = paste0("spot", 1:6)
        )
    )
    cdf <- colData(experiments(mase)[["rna"]])
    out <- subsetByColumn(mase, list(rna = cdf$tissue == "tumor" & cdf$nCount >= 150))
    expect_s4_class(out, "MultiAssaySpatialExperiment")
    expect_equal(colnames(experiments(out)[["rna"]]), c("spot2", "spot3"))
    expect_equal(nrow(spatialMap(out)), 2L)
    expect_equal(nrow(spatialPoints(out)[["tissue1"]]), 2L)
    expect_equal(spatialPoints(out)[["tissue1"]][["instance_id"]], c("spot2", "spot3"))
})

test_that("subsetByPolygon filters points", {
    pts <- DataFrame(x = c(1, 2, 3, 4, 5), y = c(1, 2, 3, 4, 5),
        instance_id = paste0("S", 1:5))
    mase <- MultiAssaySpatialExperiment(
        experiments = ExperimentList(
            assay1 = matrix(rnorm(25), 5, 5, dimnames = list(paste0("G", 1:5), paste0("S", 1:5)))
        ),
        colData = DataFrame(row.names = paste0("S", 1:5)),
        sampleMap = DataFrame(
            assay = factor("assay1"),
            primary = paste0("S", 1:5),
            colname = paste0("S", 1:5)
        ),
        points = PointsLayerList(coords = pts),
        spatialMap = DataFrame(
            assay = factor("assay1"),
            colname = paste0("S", 1:5),
            element_type = "points",
            region = "coords",
            instance_id = paste0("S", 1:5)
        )
    )
    poly <- sf::st_polygon(list(matrix(c(1.5, 1.5, 4.5, 1.5, 4.5, 4.5, 1.5, 4.5, 1.5, 1.5), ncol = 2, byrow = TRUE)))
    out <- subsetByPolygon(mase, poly)
    expect_s4_class(out, "MultiAssaySpatialExperiment")
    expect_equal(nrow(spatialPoints(out)[["coords"]]), 3L)
    expect_equal(ncol(experiments(out)[["assay1"]]), 3L)
})
