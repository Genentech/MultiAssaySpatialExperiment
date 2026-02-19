test_that("RasterLayerList constructor and validity", {
    rel <- RasterLayerList()
    expect_s4_class(rel, "RasterLayerList")
    expect_length(rel, 0L)
    expect_true(validObject(rel))

    mat <- matrix(1:6, 2, 3)
    rel2 <- RasterLayerList(img1 = mat)
    expect_length(rel2, 1L)
    expect_equal(names(rel2), "img1")
    expect_true(validObject(rel2))

    arr <- array(1:24, c(2, 3, 4))
    rel3 <- RasterLayerList(img1 = mat, img2 = arr)
    expect_length(rel3, 2L)
    expect_true(validObject(rel3))
})

test_that("RasterLayerList invalid elements fail", {
    expect_error(validObject(RasterLayerList(x = 1:5)), "dim\\(\\)")
})

test_that("RasterLayerList requires named elements", {
    mat <- matrix(1:6, 2, 3)
    expect_error(validObject(RasterLayerList(mat)), "named")
})

test_that("PointsLayerList constructor and validity", {
    pel <- PointsLayerList()
    expect_s4_class(pel, "PointsLayerList")
    expect_true(validObject(pel))

    df <- DataFrame(x = c(1, 2), y = c(3, 4))
    pel2 <- PointsLayerList(pts = df)
    expect_length(pel2, 1L)
    expect_true(validObject(pel2))
})

test_that("PointsLayerList invalid elements fail", {
    df <- DataFrame(x = 1)
    expect_error(validObject(PointsLayerList(pts = df)), "at least 2")
})

test_that("ShapesLayerList constructor and validity", {
    sel <- ShapesLayerList()
    expect_s4_class(sel, "ShapesLayerList")
    expect_true(validObject(sel))

    df <- DataFrame(id = 1:2, value = c("a", "b"))
    sel2 <- ShapesLayerList(regions = df)
    expect_length(sel2, 1L)
    expect_true(validObject(sel2))

    ## data.frame is coerced to DataFrame
    tbl <- data.frame(id = 1:2, value = c("a", "b"))
    sel3 <- ShapesLayerList(regions = tbl)
    expect_s4_class(sel3[[1L]], "DataFrame")
    expect_true(validObject(sel3))
})

test_that("Constructor identity return and guard", {
    rel <- RasterLayerList(img = matrix(1:4, 2, 2))
    expect_identical(RasterLayerList(rel), rel)

    pel <- PointsLayerList(pts = DataFrame(x = 1, y = 2))
    expect_identical(PointsLayerList(pel), pel)

    sel <- ShapesLayerList(r = DataFrame(a = 1))
    expect_identical(ShapesLayerList(sel), sel)

    expect_error(RasterLayerList(MultiAssaySpatialExperiment()), "Did you mean")
})

test_that("isEmpty and show work", {
    expect_true(isEmpty(RasterLayerList()))
    expect_true(isEmpty(PointsLayerList()))
    expect_true(isEmpty(ShapesLayerList()))

    rel <- RasterLayerList(img = matrix(1:4, 2, 2))
    expect_false(isEmpty(rel))

    expect_output(show(rel), "RasterLayerList")
    expect_output(show(rel), "img")
})
