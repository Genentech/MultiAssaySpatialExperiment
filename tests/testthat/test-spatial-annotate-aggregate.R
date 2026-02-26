### =========================================================================
### test-spatial-annotate-aggregate.R — tests for annotateWithRegions,
### aggregateByRegion
### -------------------------------------------------------------------------

test_that("annotateWithRegions adds cells column to spatialMap", {
    pts <- S4Vectors::DataFrame(
        x = c(1.5, 2.5, 3.5),
        y = c(1.5, 2.5, 3.5),
        instance_id = c("A", "B", "C"))
    shp_df <- S4Vectors::DataFrame(
        instance_id = c("cell1", "cell2", "cell3"),
        geometry = sf::st_sfc(
            sf::st_polygon(list(matrix(c(1,1,2,1,2,2,1,2,1,1), ncol=2, byrow=TRUE))),
            sf::st_polygon(list(matrix(c(2,1,3,1,3,2,2,2,2,1), ncol=2, byrow=TRUE))),
            sf::st_polygon(list(matrix(c(2,2,3,2,3,3,2,3,2,2), ncol=2, byrow=TRUE)))))
    mase <- MultiAssaySpatialExperiment(
        experiments = ExperimentList(assay1 = matrix(1:9, 3, 3,
            dimnames = list(paste0("G", 1:3), c("A", "B", "C")))),
        colData = S4Vectors::DataFrame(row.names = "s1"),
        sampleMap = S4Vectors::DataFrame(
            assay = "assay1", primary = "s1", colname = c("A", "B", "C")),
        points = PointsLayerList(centroids = pts),
        shapes = ShapesLayerList(cells = shp_df),
        spatialMap = S4Vectors::DataFrame(
            assay = "assay1", colname = c("A", "B", "C"),
            element_type = "points", region = "centroids", instance_id = c("A", "B", "C")))
    mase <- annotateWithRegions(mase, points = "centroids", shapes = "cells")
    spmap <- spatialMap(mase)
    expect_true("cells" %in% colnames(spmap))
    expect_equal(spmap[["cells"]][1L], "cell1")
    expect_equal(spmap[["cells"]][2L], "cell3")
    expect_true(is.na(spmap[["cells"]][3L]))
})

test_that("aggregateByRegion sums assay values by shape", {
    pts <- S4Vectors::DataFrame(
        x = c(1.5, 2.5, 3.5),
        y = c(1.5, 2.5, 3.5),
        instance_id = c("A", "B", "C"))
    shp_df <- S4Vectors::DataFrame(
        instance_id = c("cell1", "cell2", "cell3"),
        geometry = sf::st_sfc(
            sf::st_polygon(list(matrix(c(1,1,2,1,2,2,1,2,1,1), ncol=2, byrow=TRUE))),
            sf::st_polygon(list(matrix(c(2,1,3,1,3,2,2,2,2,1), ncol=2, byrow=TRUE))),
            sf::st_polygon(list(matrix(c(2,2,3,2,3,3,2,3,2,2), ncol=2, byrow=TRUE)))))
    mase <- MultiAssaySpatialExperiment(
        experiments = ExperimentList(assay1 = matrix(1:9, 3, 3,
            dimnames = list(paste0("G", 1:3), c("A", "B", "C")))),
        colData = S4Vectors::DataFrame(row.names = "s1"),
        sampleMap = S4Vectors::DataFrame(
            assay = "assay1", primary = "s1", colname = c("A", "B", "C")),
        points = PointsLayerList(centroids = pts),
        shapes = ShapesLayerList(cells = shp_df),
        spatialMap = S4Vectors::DataFrame(
            assay = "assay1", colname = c("A", "B", "C"),
            element_type = "points", region = "centroids", instance_id = c("A", "B", "C")))
    mase <- annotateWithRegions(mase, points = "centroids", shapes = "cells")
    agg <- aggregateByRegion(mase, by = "cells", FUN = "sum")
    expect_type(agg, "list")
    expect_equal(names(agg), "assay1")
    expect_equal(agg[["assay1"]]["G1", "cell1"], 1)
    expect_equal(agg[["assay1"]]["G1", "cell3"], 4)
    expect_equal(agg[["assay1"]]["G3", "cell1"], 3)
    expect_equal(agg[["assay1"]]["G3", "cell3"], 6)
})

test_that("aggregateByRegion count returns DataFrame", {
    pts <- S4Vectors::DataFrame(
        x = c(1.5, 2.5), y = c(1.5, 2.5), instance_id = c("A", "B"))
    shp_df <- S4Vectors::DataFrame(
        instance_id = "cell1",
        geometry = sf::st_sfc(sf::st_polygon(list(matrix(c(1,1,3,1,3,3,1,3,1,1), ncol=2, byrow=TRUE)))))
    mase <- MultiAssaySpatialExperiment(
        experiments = ExperimentList(assay1 = matrix(1:6, 2, 3, dimnames = list(c("G1", "G2"), c("A", "B", "C")))),
        colData = S4Vectors::DataFrame(row.names = "s1"),
        sampleMap = S4Vectors::DataFrame(assay = "assay1", primary = "s1", colname = c("A", "B", "C")),
        points = PointsLayerList(centroids = pts),
        shapes = ShapesLayerList(cells = shp_df),
        spatialMap = S4Vectors::DataFrame(
            assay = "assay1", colname = c("A", "B"),
            element_type = "points", region = "centroids", instance_id = c("A", "B")))
    mase <- annotateWithRegions(mase, points = "centroids", shapes = "cells")
    cnt <- aggregateByRegion(mase, by = "cells", FUN = "count")
    expect_s4_class(cnt, "DataFrame")
    expect_equal(cnt$region, "cell1")
    expect_equal(cnt$count, 2L)
})

test_that("aggregateByRegion mean aggregates correctly", {
    pts <- S4Vectors::DataFrame(
        x = c(1.5, 2.5), y = c(1.5, 2.5), instance_id = c("A", "B"))
    shp_df <- S4Vectors::DataFrame(
        instance_id = "cell1",
        geometry = sf::st_sfc(sf::st_polygon(list(matrix(c(1,1,3,1,3,3,1,3,1,1), ncol=2, byrow=TRUE)))))
    mat <- matrix(c(2, 4, 6, 8), 2, 2, dimnames = list(c("G1", "G2"), c("A", "B")))
    mase <- MultiAssaySpatialExperiment(
        experiments = ExperimentList(assay1 = mat),
        colData = S4Vectors::DataFrame(row.names = "s1"),
        sampleMap = S4Vectors::DataFrame(assay = "assay1", primary = "s1", colname = c("A", "B")),
        points = PointsLayerList(centroids = pts),
        shapes = ShapesLayerList(cells = shp_df),
        spatialMap = S4Vectors::DataFrame(
            assay = "assay1", colname = c("A", "B"),
            element_type = "points", region = "centroids", instance_id = c("A", "B")))
    mase <- annotateWithRegions(mase, points = "centroids", shapes = "cells")
    agg <- aggregateByRegion(mase, by = "cells", FUN = "mean")
    expect_equal(agg[["assay1"]]["G1", "cell1"], 4)
    expect_equal(agg[["assay1"]]["G2", "cell1"], 6)
})
