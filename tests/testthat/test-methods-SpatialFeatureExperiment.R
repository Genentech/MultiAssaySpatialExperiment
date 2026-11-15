### =========================================================================
### Test methods for SpatialFeatureExperiment generics (Phase 7)
### =========================================================================

test_that("SFE methods work via ephemeral reconstruction from a plain SPE-derived MASE", {
    skip_if_not_installed("SpatialFeatureExperiment")
    library(SpatialFeatureExperiment)
    mat <- matrix(1:12, 3, 4)
    colnames(mat) <- paste0("cell", 1:4)
    spe <- SpatialExperiment::SpatialExperiment(
        list(counts = mat),
        spatialCoords = matrix(rnorm(8), 4, 2)
    )
    mase <- as(spe, "MultiAssaySpatialExperiment")
    u <- unit(mase)
    expect_true(is.character(u) || is.null(u))
})

test_that("SFE methods delegate when SFE assay present", {
    skip_if_not_installed("SpatialFeatureExperiment")
    library(SpatialFeatureExperiment)
    mat <- matrix(1:12, 3, 4)
    colnames(mat) <- paste0("cell", 1:4)
    spe <- SpatialExperiment::SpatialExperiment(
        list(counts = mat),
        spatialCoords = matrix(rnorm(8), 4, 2)
    )
    sfe <- as(spe, "SpatialFeatureExperiment")
    mase <- as(sfe, "MultiAssaySpatialExperiment")
    u <- unit(mase)
    expect_true(is.character(u) || is.null(u))
})

test_that("dimGeometry errors when the MASE has more than one assay", {
    skip_if_not_installed("SpatialFeatureExperiment")
    library(SpatialFeatureExperiment)
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
    expect_error(dimGeometry(mase, MARGIN = 2),
        "requires exactly one assay"
    )
})

test_that("bbox delegates to SFE when present", {
    skip_if_not_installed("SpatialFeatureExperiment")
    library(SpatialFeatureExperiment)
    mat <- matrix(1:12, 3, 4)
    colnames(mat) <- paste0("cell", 1:4)
    spe <- SpatialExperiment::SpatialExperiment(
        list(counts = mat),
        spatialCoords = matrix(rnorm(8), 4, 2, dimnames = list(NULL, c("x", "y")))
    )
    sfe <- as(spe, "SpatialFeatureExperiment")
    mase <- as(sfe, "MultiAssaySpatialExperiment")
    b <- bbox(mase)
    expect_true(is.numeric(b))
    expect_true(length(b) >= 4L)
})

.make_populated_sfe <- function() {
    mat <- matrix(rpois(24, 3), 4, 6,
        dimnames = list(paste0("g", 1:4), paste0("cell", 1:6)))
    coords <- matrix(c(0, 0, 1, 0, 2, 0, 0, 1, 1, 1, 2, 1), 6, 2,
        byrow = TRUE, dimnames = list(NULL, c("x", "y")))
    spe <- SpatialExperiment::SpatialExperiment(list(counts = mat),
        spatialCoords = coords)
    sfe <- as(spe, "SpatialFeatureExperiment")

    # a second, non-centroid colGeometry: a small square around each cell
    cg2 <- sf::st_sf(geometry = sf::st_sfc(lapply(seq_len(6), function(i) {
        cx <- coords[i, 1]; cy <- coords[i, 2]
        sf::st_polygon(list(matrix(c(
            cx - .4, cy - .4, cx + .4, cy - .4,
            cx + .4, cy + .4, cx - .4, cy + .4,
            cx - .4, cy - .4), 5, 2, byrow = TRUE)))
    })))
    rownames(cg2) <- colnames(sfe)
    SpatialFeatureExperiment::colGeometry(sfe, "footprint") <- cg2

    # an annotGeometry: one region covering every cell (not one row per cell)
    ag <- sf::st_sf(sample_id = "sample01",
        geometry = sf::st_sfc(sf::st_polygon(list(matrix(c(
            -1, -1, 3, -1, 3, 2, -1, 2, -1, -1), 5, 2, byrow = TRUE)))))
    SpatialFeatureExperiment::annotGeometry(sfe, "tissue") <- ag

    # a col-margin spatial neighbor graph
    g <- SpatialFeatureExperiment::findSpatialNeighbors(sfe, MARGIN = 2L,
        method = "knearneigh", k = 2)
    SpatialFeatureExperiment::colGraph(sfe, "knn") <- g
    sfe
}

test_that("SFE -> MASE -> SFE round-trips multiple colGeometries and an annotGeometry", {
    skip_if_not_installed("SpatialFeatureExperiment")
    library(SpatialFeatureExperiment)
    sfe <- .make_populated_sfe()

    mase <- as(sfe, "MultiAssaySpatialExperiment")
    expect_true(validObject(mase))
    expect_setequal(spatialShapeNames(mase), c("footprint", "tissue"))
    expect_setequal(as.character(spatialMap(mase)[["region"]]),
        c("coordinates", "footprint"))

    sfe2 <- as(mase, "SpatialFeatureExperiment")
    expect_true(validObject(sfe2))
    expect_setequal(colGeometryNames(sfe2), c("centroids", "footprint"))
    expect_setequal(annotGeometryNames(sfe2), "tissue")
    expect_equal(ncol(sfe2), ncol(sfe))
})

test_that("an SFE carrying the legacy 'full_res_image_pixels' unit round-trips through MASE", {
    skip_if_not_installed("SpatialFeatureExperiment")
    library(SpatialFeatureExperiment)
    sfe <- .make_populated_sfe()
    SingleCellExperiment::int_metadata(sfe)$unit <- "full_res_image_pixels"
    expect_identical(unit(sfe), "full_res_image_pixels")

    mase <- as(sfe, "MultiAssaySpatialExperiment")
    expect_identical(S4Vectors::metadata(mase)[["spatialUnit"]], "full_res_image_pixel")

    sfe2 <- as(mase, "SpatialFeatureExperiment")
    expect_true(validObject(sfe2))
    expect_identical(unit(sfe2), "full_res_image_pixel")
    expect_equal(ncol(sfe2), ncol(sfe))

    # a MASE written before this fix still carries the legacy spelling
    S4Vectors::metadata(mase)[["spatialUnit"]] <- "full_res_image_pixels"
    expect_identical(unit(as(mase, "SpatialFeatureExperiment")), "full_res_image_pixel")
})

test_that("col-margin spatialGraphs round-trip through MASE natively and via SFE", {
    skip_if_not_installed("SpatialFeatureExperiment")
    library(SpatialFeatureExperiment)
    sfe <- .make_populated_sfe()
    g <- colGraph(sfe, type = "knn")

    mase <- as(sfe, "MultiAssaySpatialExperiment")
    expect_setequal(spatialGraphNames(mase, MARGIN = 2L), "knn")
    g2 <- spatialGraphs(mase, MARGIN = 2L, name = "knn")
    expect_s3_class(g2, "listw")
    expect_equal(unname(as.list(g2[["neighbours"]])),
        unname(as.list(g[["neighbours"]])))

    sfe2 <- as(mase, "SpatialFeatureExperiment")
    expect_setequal(colGraphNames(sfe2), "knn")
    g3 <- colGraph(sfe2, type = "knn")
    expect_equal(unname(as.list(g3[["neighbours"]])),
        unname(as.list(g[["neighbours"]])))
})

test_that("dimGeometry/aggregate/findDebrisCells work against a real populated SFE-derived MASE", {
    skip_if_not_installed("SpatialFeatureExperiment")
    library(SpatialFeatureExperiment)
    sfe <- .make_populated_sfe()
    mase <- as(sfe, "MultiAssaySpatialExperiment")

    fp <- dimGeometry(mase, type = "footprint", MARGIN = 2L)
    expect_equal(nrow(fp), 6L)

    region_sf <- annotGeometry(sfe, "tissue")
    agg <- aggregate(mase, by = region_sf)
    expect_s4_class(agg, "MultiAssaySpatialExperiment")

    dc <- findDebrisCells(mase)
    expect_false(is.null(dc))
})
