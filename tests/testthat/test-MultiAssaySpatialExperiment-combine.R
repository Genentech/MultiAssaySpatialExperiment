library(MultiAssayExperiment)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### c()
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.simpleMASE <- function(assay_nm = "assay1", primary_prefix = "S", n = 4L) {
    mat <- matrix(rnorm(5L * n), nrow = 5L, ncol = n,
        dimnames = list(paste0("G", 1:5), paste0(primary_prefix, seq_len(n))))
    cd <- DataFrame(row.names = paste0(primary_prefix, seq_len(n)))
    sm <- DataFrame(
        assay = factor(assay_nm, assay_nm),
        primary = paste0(primary_prefix, seq_len(n)),
        colname = paste0(primary_prefix, seq_len(n))
    )
    MultiAssaySpatialExperiment(
        experiments = ExperimentList(setNames(list(mat), assay_nm)),
        colData = cd,
        sampleMap = sm
    )
}

.simpleMASE_with_spatial <- function(assay_nm = "assay1", layer_nm = "pts",
    primary_prefix = "S", n = 4L) {
    mase <- .simpleMASE(assay_nm = assay_nm, primary_prefix = primary_prefix, n = n)
    pts <- DataFrame(x = seq_len(n), y = seq_len(n),
        instance_id = paste0(primary_prefix, seq_len(n)))
    mase@points <- PointsLayerList(setNames(list(pts), layer_nm))
    mase
}

test_that("c() merges two MultiAssaySpatialExperiment objects", {
    mase1 <- .simpleMASE("assay1", "S", 4L)
    mase2 <- .simpleMASE("assay2", "S", 4L)  # same primaries, different assay
    combined <- c(mase1, mase2)
    expect_true(validObject(combined))
    expect_equal(names(combined), c("assay1", "assay2"))
    expect_equal(nrow(colData(combined)), 4L)
    expect_equal(nrow(sampleMap(combined)),
        nrow(sampleMap(mase1)) + nrow(sampleMap(mase2)))
})

test_that("c() requires unique experiment names", {
    mase1 <- .simpleMASE("assay1")
    mase2 <- .simpleMASE("assay1")
    expect_error(c(mase1, mase2), "unique experiment names")
})

test_that("c() requires unique spatial layer names when both have layers", {
    mase1 <- .simpleMASE_with_spatial("assay1", "pts1")
    mase2 <- .simpleMASE_with_spatial("assay2", "pts1")
    expect_error(c(mase1, mase2), "unique spatial layer names")
})

test_that("c() merges spatial layers with unique names", {
    mase1 <- .simpleMASE_with_spatial("assay1", "pts_x")
    mase2 <- .simpleMASE_with_spatial("assay2", "pts_y")
    combined <- c(mase1, mase2)
    expect_equal(names(spatialPoints(combined)), c("pts_x", "pts_y"))
    expect_true(validObject(combined))
})

test_that("c() delegates to MAE when adding experiments (not MASE)", {
    mase <- .simpleMASE("assay1")
    mat2 <- matrix(rnorm(20), nrow = 5, ncol = 4,
        dimnames = list(paste0("G", 1:5), paste0("S", 1:4)))
    result <- c(mase, assay2 = mat2, mapFrom = 1L)
    expect_true(validObject(result))
    expect_equal(names(result), c("assay1", "assay2"))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### complete.cases()
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

test_that("complete.cases() identifies samples in all assays", {
    mat1 <- matrix(rnorm(20), nrow = 5, ncol = 4,
        dimnames = list(paste0("G", 1:5), paste0("S", 1:4)))
    mat2 <- matrix(rnorm(15), nrow = 5, ncol = 3,
        dimnames = list(paste0("G", 1:5), paste0("S", 1:3)))
    cd <- DataFrame(row.names = paste0("S", 1:4))
    sm <- rbind(
        DataFrame(assay = "assay1", primary = paste0("S", 1:4), colname = paste0("S", 1:4)),
        DataFrame(assay = "assay2", primary = paste0("S", 1:3), colname = paste0("S", 1:3)))
    mase <- MultiAssaySpatialExperiment(
        experiments = ExperimentList(assay1 = mat1, assay2 = mat2),
        colData = cd,
        sampleMap = sm
    )
    cc <- complete.cases(mase)
    expect_true(is.logical(cc))
    expect_equal(length(cc), 4L)
    expect_equal(sum(cc), 3L)
    expect_false(cc[[4L]])
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### cbind()
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

test_that("cbind() returns single object unchanged", {
    mase <- .simpleMASE()
    expect_identical(cbind(mase), mase)
})

test_that("cbind() combines two MASE with same assay names", {
    mase1 <- .simpleMASE(primary_prefix = "A", n = 2L)
    mase2 <- .simpleMASE(primary_prefix = "B", n = 2L)
    combined <- cbind(mase1, mase2)
    expect_true(validObject(combined))
    expect_equal(ncol(combined[["assay1"]]), 4L)
    expect_equal(nrow(colData(combined)), 4L)
    expect_equal(rownames(colData(combined)), c("A1", "A2", "B1", "B2"))
})

test_that("cbind() requires same assay names", {
    mase1 <- .simpleMASE("assay1")
    mase2 <- .simpleMASE("assay2")
    expect_error(cbind(mase1, mase2), "same assay names")
})

test_that("cbind() merges spatial point layers with same name", {
    mase1 <- .simpleMASE_with_spatial(layer_nm = "centroids", primary_prefix = "A", n = 2L)
    mase2 <- .simpleMASE_with_spatial(layer_nm = "centroids", primary_prefix = "B", n = 2L)
    combined <- cbind(mase1, mase2)
    expect_equal(nrow(spatialPoints(combined)[["centroids"]]), 4L)
    expect_true(validObject(combined))
})
