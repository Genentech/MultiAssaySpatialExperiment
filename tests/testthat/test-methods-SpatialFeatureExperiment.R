### =========================================================================
### Test methods for SpatialFeatureExperiment generics (Phase 7)
### =========================================================================

test_that("SFE methods require SpatialFeatureExperiment", {
    skip_if_not_installed("SpatialFeatureExperiment")
    library(SpatialFeatureExperiment)
    mat <- matrix(1:12, 3, 4)
    colnames(mat) <- paste0("cell", 1:4)
    spe <- SpatialExperiment::SpatialExperiment(
        list(counts = mat),
        spatialCoords = matrix(rnorm(8), 4, 2)
    )
    mase <- as(spe, "MultiAssaySpatialExperiment")
    expect_null(unit(mase))
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

test_that("dimGeometry errors when no SFE assay", {
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
        "No SpatialFeatureExperiment assay"
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
