test_that("buildSpatialMap assembles FK rows from sampleMap", {
    sm <- DataFrame(
        assay = factor(rep("rna", 3), "rna"),
        primary = c("S1", "S1", "S2"),
        colname = paste0("spot", 1:3)
    )
    spmap <- buildSpatialMap(sm, region = "tissue1", element_type = "points")
    expect_s4_class(spmap, "DataFrame")
    expect_equal(ncol(spmap), 5L)
    expect_equal(nrow(spmap), 3L)
    expect_true(all(spmap[["region"]] == "tissue1"))
    expect_equal(spmap[["instance_id"]], paste0("spot", 1:3))
})

test_that("prepMASE harmonizes spatialMap with sampleMap", {
    mat <- matrix(1:6, 2, 3, dimnames = list(c("G1", "G2"), paste0("spot", 1:3)))
    sm <- DataFrame(
        assay = factor(rep("rna", 3), "rna"),
        primary = c("S1", "S1", "S2"),
        colname = paste0("spot", 1:3)
    )
    pts <- DataFrame(x = 1:3, y = 1:3, instance_id = paste0("spot", 1:3))
    spmap <- buildSpatialMap(sm, "tissue1", "points")
    orphan <- DataFrame(
        assay = factor("rna", "rna"),
        colname = "missing",
        element_type = "points",
        region = "tissue1",
        instance_id = "ghost"
    )
    prepared <- prepMASE(
        ExperimentList(rna = mat),
        DataFrame(row.names = c("S1", "S2")),
        sm,
        points = PointsLayerList(tissue1 = pts),
        spatialMap = rbind(spmap, orphan)
    )
    expect_equal(nrow(prepared$spatialMap), 3L)
    expect_true(length(prepared$metadata$drops) > 0L)
    mase <- do.call(MultiAssaySpatialExperiment, prepared)
    expect_s4_class(mase, "MultiAssaySpatialExperiment")
})

test_that("checkSpatialMap reports FK problems", {
    sm <- DataFrame(
        assay = factor("rna", "rna"),
        primary = "S1",
        colname = "spot1"
    )
    spmap <- buildSpatialMap(sm, "coords", "points")
    spmap[["instance_id"]] <- "missing"
    msgs <- checkSpatialMap(spmap, sm, points = PointsLayerList(
        coords = DataFrame(x = 1, y = 1, instance_id = "spot1")
    ))
    expect_true(length(msgs) > 0L)
})
