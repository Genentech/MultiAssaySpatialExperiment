test_that("labelsToShapes and shapesToLabels round-trip label values", {
    skip_if_not_installed("terra")
    skip_if_not_installed("sf")

    mat <- matrix(c(0L, 1L, 1L, 2L), 2, 2)
    polys <- labelsToShapes(mat, dissolve = TRUE)
    expect_s3_class(polys, "sf")
    expect_gt(nrow(polys), 0L)

    polys$instance_id <- polys[["lyr.1"]]
    labels <- shapesToLabels(polys, field = "instance_id", background = 0L)
    expect_s4_class(labels, "SpatRaster")
})

test_that("readers attach loaded image rasters to spatialImages", {
    skip_if_not_installed("DropletUtils")
    skip_if_not_installed("Matrix")
    skip_if_not_installed("jsonlite")
    skip_if_not_installed("sf")
    skip_if_not_installed("terra")

    data_dir <- make_mock_visium_dir(n_spots = 2L)
    on.exit(unlink(data_dir, recursive = TRUE), add = TRUE)

    spatial_dir <- file.path(data_dir, "spatial")
    png_files <- list.files(spatial_dir, pattern = "\\.png$", full.names = TRUE)
    for (old_png in png_files)
        file.remove(old_png)

    png_path <- file.path(spatial_dir, "tissue_lowres_image.png")
    terra::writeRaster(
        terra::rast(matrix(1:12, 3, 4)),
        png_path,
        overwrite = TRUE
    )

    mase <- readVisiumMASE(data_dir, images = TRUE, load_images = TRUE)
    expect_s4_class(mase, "MultiAssaySpatialExperiment")
    expect_gt(length(spatialImages(mase)), 0L)
    expect_true("tissue_lowres_image" %in% spatialImageNames(mase))
})

test_that(".spatialImages_from_imgData skips NULL raster payloads", {
    imgData <- DataFrame(
        sample_id = "s1",
        image_id = "img1",
        data = I(list(NULL)),
        scaleFactor = 1.0,
        width = 10L,
        height = 10L,
        path = "dummy.png"
    )
    imgs <- MultiAssaySpatialExperiment:::.spatialImages_from_imgData(imgData)
    expect_s4_class(imgs, "RasterLayerList")
    expect_equal(length(imgs), 0L)
})
