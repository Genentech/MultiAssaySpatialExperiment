## Register SpatialData coercion when the optional package loads.
.has_spatialdata_coercion <- function() {
    if (!requireNamespace("spatialdataR", quietly = TRUE))
        return(FALSE)
    if (!isNamespaceLoaded("spatialdataR"))
        loadNamespace("spatialdataR")
    tryCatch({
        methods::getMethod("coerce", c("SpatialData", "MultiAssaySpatialExperiment"))
        TRUE
    }, error = function(e) FALSE)
}

## Unit squares as an sf: spatialdataR shapes need real geometry.
.test_polys <- function(...) {
    d <- data.frame(...)
    sq <- lapply(seq_len(nrow(d)), function(i)
        sf::st_polygon(list(rbind(c(i, 0), c(i + 1, 0), c(i + 1, 1),
                                  c(i, 1), c(i, 0)))))
    sf::st_sf(d, geometry = sf::st_sfc(sq))
}

test_that("as(SpatialData, 'MultiAssaySpatialExperiment') works when SpatialData available", {
    skip_if_not_installed("spatialdataR")
    expect_true(.has_spatialdata_coercion(),
                info = "SpatialData coercion not registered")
    library(SingleCellExperiment)
    sce <- SingleCellExperiment(
        list(counts = matrix(1:12, 3, 4)),
        colData = DataFrame(row.names = paste0("cell", 1:4))
    )
    sd <- spatialdataR::SpatialData(
        images = list(),
        labels = list(),
        points = list(
            pts1 = spatialdataR::SpatialDataPoint(data = data.frame(
                x = c(1, 2, 3, 4),
                y = c(1, 2, 3, 4),
                instance_id = 1:4
            ))
        ),
        shapes = list(
            shp1 = spatialdataR::SpatialDataShape(data = .test_polys(
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
    skip_if_not_installed("spatialdataR")
    expect_true(.has_spatialdata_coercion(),
                info = "SpatialData coercion not registered")
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
    expect_equal(length(spatialdataR::tables(sd)), 1L)
    expect_equal(spatialdataR::tableNames(sd), "assay1")
    expect_equal(spatialdataR::pointNames(sd), "coords")
})

test_that("SD-MASE-SD round-trip preserves tables and points", {
    skip_if_not_installed("spatialdataR")
    expect_true(.has_spatialdata_coercion(),
                info = "SpatialData coercion not registered")
    library(SingleCellExperiment)
    sce <- SingleCellExperiment(
        list(counts = matrix(1:12, 3, 4)),
        colData = DataFrame(row.names = paste0("cell", 1:4))
    )
    sd1 <- spatialdataR::SpatialData(
        images = list(),
        labels = list(),
        points = list(
            p1 = spatialdataR::SpatialDataPoint(data = data.frame(
                x = 1:4, y = 2:5, id = 1:4
            ))
        ),
        shapes = list(),
        tables = list(t1 = sce)
    )
    mase <- as(sd1, "MultiAssaySpatialExperiment")
    sd2 <- as(mase, "SpatialData")
    expect_equal(spatialdataR::tableNames(sd2), "t1")
    expect_equal(ncol(spatialdataR::table(sd2, "t1")), 4L)
})

test_that("multi-table SpatialData coerces to multi-assay MASE", {
    skip_if_not_installed("spatialdataR")
    expect_true(.has_spatialdata_coercion(),
                info = "SpatialData coercion not registered")
    library(SingleCellExperiment)
    sce1 <- SingleCellExperiment(
        list(counts = matrix(1:6, 2, 3)),
        colData = DataFrame(row.names = paste0("A", 1:3))
    )
    sce2 <- SingleCellExperiment(
        list(counts = matrix(7:12, 2, 3)),
        colData = DataFrame(row.names = paste0("B", 1:3))
    )
    sd <- spatialdataR::SpatialData(
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
    skip_if_not_installed("spatialdataR")
    expect_true(.has_spatialdata_coercion(),
                info = "SpatialData coercion not registered")
    sd <- spatialdataR::SpatialData(
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
    skip_if_not_installed("spatialdataR")
    expect_true(.has_spatialdata_coercion(),
                info = "SpatialData coercion not registered")
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
    sd <- spatialdataR::SpatialData(
        images = list(),
        labels = list(),
        points = list(
            cells = spatialdataR::SpatialDataPoint(data = data.frame(
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
    skip_if_not_installed("spatialdataR")
    expect_true(.has_spatialdata_coercion(),
                info = "SpatialData coercion not registered")
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
    sd <- spatialdataR::SpatialData(
        images = list(),
        labels = list(),
        points = list(),
        shapes = list(
            nuclei = spatialdataR::SpatialDataShape(data = .test_polys(
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
    skip_if_not_installed("spatialdataR")
    expect_true(.has_spatialdata_coercion(),
                info = "SpatialData coercion not registered")
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

.blobs_zarr <- function() system.file("extdata", "blobs.zarr", package = "spatialdataR")

test_that("a real SpatialData converts its tables, points and shapes and reports the dropped images and labels", {
    skip_if_not_installed("spatialdataR")
    expect_true(.has_spatialdata_coercion(),
                info = "SpatialData coercion not registered")
    skip_if(!nzchar(.blobs_zarr()), "spatialdataR example data not found")
    sd <- spatialdataR::readSpatialData(.blobs_zarr())

    expect_warning(
        mase <- as(sd, "MultiAssaySpatialExperiment"),
        "images and labels are not converted.*blobs_image.*blobs_labels")
    expect_setequal(spatialPointNames(mase), spatialdataR::pointNames(sd))
    expect_setequal(spatialShapeNames(mase), spatialdataR::shapeNames(sd))
    # the lazy points are materialized, not carried as a query
    expect_equal(nrow(spatialPoints(mase)[["blobs_points"]]), 200L)
    expect_length(spatialImages(mase), 0L)
    expect_length(spatialLabels(mase), 0L)
})

test_that("a table that annotates a dropped label loses that annotation on the way back", {
    skip_if_not_installed("spatialdataR")
    expect_true(.has_spatialdata_coercion(),
                info = "SpatialData coercion not registered")
    skip_if(!nzchar(.blobs_zarr()), "spatialdataR example data not found")
    sd <- spatialdataR::readSpatialData(.blobs_zarr())
    region <- SingleCellExperiment::int_metadata(
        spatialdataR::tables(sd)[[1L]])$spatialdata_attrs$region
    expect_true(region %in% spatialdataR::labelNames(sd))

    mase <- suppressWarnings(as(sd, "MultiAssaySpatialExperiment"))
    expect_warning(back <- as(mase, "SpatialData"),
                   paste0("annotates ", region, ".*not converted"))
    expect_s4_class(back, "SpatialData")
    expect_setequal(spatialdataR::pointNames(back), spatialdataR::pointNames(sd))
    expect_setequal(spatialdataR::shapeNames(back), spatialdataR::shapeNames(sd))
    expect_null(SingleCellExperiment::int_metadata(
        spatialdataR::tables(back)[[1L]])$spatialdata_attrs)
})

test_that("MASE images and labels are reported as dropped when converting to SpatialData", {
    skip_if_not_installed("spatialdataR")
    expect_true(.has_spatialdata_coercion(),
                info = "SpatialData coercion not registered")
    library(SingleCellExperiment)
    sce <- SingleCellExperiment(
        list(counts = matrix(1:12, 3, 4)),
        colData = DataFrame(row.names = paste0("cell", 1:4)))
    mase <- MultiAssaySpatialExperiment(
        experiments = ExperimentList(assay1 = sce),
        colData = DataFrame(row.names = "sample01"),
        sampleMap = DataFrame(
            assay = factor("assay1", "assay1"),
            primary = "sample01",
            colname = paste0("cell", 1:4)),
        images = RasterLayerList(brightfield = matrix(1:12, 3, 4)))
    expect_warning(sd <- as(mase, "SpatialData"),
                   "images and labels are not converted.*brightfield")
    expect_s4_class(sd, "SpatialData")
})

test_that("the integration works in a clean session, whichever package loads first", {
    skip_if_not_installed("spatialdataR")
    rscript <- file.path(R.home("bin"), "Rscript")
    body <- '
        mase <- MultiAssaySpatialExperiment::MultiAssaySpatialExperiment(
            experiments = MultiAssayExperiment::ExperimentList(
                a = SingleCellExperiment::SingleCellExperiment(list(
                    counts = matrix(1:4, 2, dimnames = list(c("g1", "g2"), c("s1", "s2")))))),
            colData = S4Vectors::DataFrame(row.names = c("s1", "s2")),
            sampleMap = S4Vectors::DataFrame(
                assay = factor(c("a", "a"), "a"), primary = c("s1", "s2"),
                colname = c("s1", "s2")),
            points = MultiAssaySpatialExperiment::PointsLayerList(
                p = S4Vectors::DataFrame(x = 1:2, y = 1:2, instance_id = c("s1", "s2"))))
        sd <- as(mase, "SpatialData")
        back <- as(sd, "MultiAssaySpatialExperiment")
        stopifnot(identical(spatialdataR::layer(mase, "p"), "points"),
                  identical(spatialdataR::element(mase, 1L),
                            spatialdataR::element(mase, "p")),
                  identical(spatialdataR::pointNames(sd), "p"))
        # a table that annotates a points element must be linked to it
        sce <- SingleCellExperiment::SingleCellExperiment(
            list(counts = matrix(1:12, 3, 4)),
            colData = S4Vectors::DataFrame(
                rk = rep("cells", 4), instance_id = paste0("cell", 1:4),
                row.names = paste0("cell", 1:4)))
        imd <- SingleCellExperiment::int_metadata(sce)
        imd[["spatialdata_attrs"]] <- list(
            region = "cells", region_key = "rk", instance_key = "instance_id")
        sce <- SingleCellExperiment::`int_metadata<-`(sce, value = imd)
        sd2 <- spatialdataR::SpatialData(
            points = list(cells = spatialdataR::SpatialDataPoint(data = data.frame(
                x = 1:4, y = 2:5, instance_id = paste0("cell", 1:4)))),
            tables = list(t1 = sce))
        sm <- MultiAssaySpatialExperiment::spatialMap(
            as(sd2, "MultiAssaySpatialExperiment"))
        stopifnot(!is.null(sm), nrow(sm) == 4L)
        cat("CLEAN_SESSION_OK")'
    # In a clean session neither SingleCellExperiment nor (second case)
    # spatialdataR is attached, which the other tests always have.
    for (setup in c('library(spatialdataR); library(MultiAssaySpatialExperiment)',
                    'library(MultiAssaySpatialExperiment); loadNamespace("spatialdataR")')) {
        code <- paste0("suppressMessages({", setup, "})\n", body)
        out <- suppressWarnings(system2(rscript, c("-e", shQuote(code)),
                                        stdout = TRUE, stderr = TRUE))
        expect_match(paste(out, collapse = "\n"), "CLEAN_SESSION_OK", fixed = TRUE,
                     info = setup)
    }
})


## A small MASE over two cells-by-table with the given layers.
.mase_with <- function(..., spatialMap = NULL, cells = c("a", "b", "c", "d")) {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(seq_along(cells) * 2L, 2, length(cells),
                             dimnames = list(c("g1", "g2"), cells))))
    MultiAssaySpatialExperiment(
        experiments = ExperimentList(t1 = sce),
        colData = DataFrame(row.names = "s"),
        sampleMap = DataFrame(assay = factor(rep("t1", length(cells)), "t1"),
                              primary = "s", colname = cells),
        ..., spatialMap = spatialMap)
}

.shapes_df <- function(n, offset = 0) {
    df <- DataFrame(id = paste0("sh", seq_len(n)))
    df$geometry <- sf::st_sfc(lapply(seq_len(n), function(i)
        sf::st_point(c(i + offset, i * 10 + offset))))
    df
}

test_that("every points and shapes element keeps its own rows and values", {
    skip_if_not_installed("spatialdataR")
    expect_true(.has_spatialdata_coercion(),
                info = "SpatialData coercion not registered")
    # spatialdataR's constructors used to share one table per kind, so each
    # element silently took on the data of the last one built
    mase <- .mase_with(
        shapes = ShapesLayerList(two = .shapes_df(2), five = .shapes_df(5, 100)),
        points = PointsLayerList(
            p3 = DataFrame(x = c(1, 2, 3), y = c(4, 5, 6), instance_id = letters[1:3]),
            p4 = DataFrame(x = c(10, 20, 30, 40), y = c(1, 2, 3, 4),
                           instance_id = letters[1:4])))
    sd <- as(mase, "SpatialData")

    sh <- function(nm) as.data.frame(spatialdataR::shape(sd, nm))
    expect_equal(nrow(sh("two")), 2L)
    expect_equal(nrow(sh("five")), 5L)
    expect_equal(sh("two")$id, paste0("sh", 1:2))
    expect_equal(sh("five")$id, paste0("sh", 1:5))
    expect_equal(unname(sf::st_coordinates(sf::st_as_sf(sh("five")))[, "X"]), 101:105)

    pt <- function(nm) as.data.frame(spatialdataR::point(sd, nm))
    expect_equal(nrow(pt("p3")), 3L)
    expect_equal(nrow(pt("p4")), 4L)
    expect_equal(pt("p3")$instance_id, letters[1:3])
    expect_equal(unname(sf::st_coordinates(sf::st_as_sf(pt("p4")))[, "X"]),
                 c(10, 20, 30, 40))
})

test_that("a SpatialData round trip keeps the rows and columns of every element", {
    skip_if_not_installed("spatialdataR")
    skip_if(!nzchar(.blobs_zarr()), "spatialdataR example data not found")
    expect_true(.has_spatialdata_coercion(),
                info = "SpatialData coercion not registered")
    sd <- spatialdataR::readSpatialData(.blobs_zarr())
    mase <- suppressWarnings(as(sd, "MultiAssaySpatialExperiment"))
    back <- suppressWarnings(as(mase, "SpatialData"))
    for (nm in c(spatialdataR::pointNames(sd), spatialdataR::shapeNames(sd))) {
        was <- as.data.frame(spatialdataR::element(sd, nm))
        now <- as.data.frame(spatialdataR::element(back, nm))
        expect_equal(nrow(now), nrow(was), info = nm)
        expect_true(all(names(was) %in% names(now)), info = nm)
    }
    # names are kept as stored, not made syntactic
    expect_true("__null_dask_index__" %in%
                    colnames(spatialPoints(mase)[["blobs_points"]]))
})

test_that("points whose coordinates are not called x and y convert", {
    skip_if_not_installed("spatialdataR")
    expect_true(.has_spatialdata_coercion(),
                info = "SpatialData coercion not registered")
    mase <- .mase_with(points = PointsLayerList(coordinates = DataFrame(
        imageX = c(5, 6, 7, 8), imageY = c(1, 2, 3, 4),
        instance_id = c("a", "b", "c", "d"))))
    sd <- as(mase, "SpatialData")
    xy <- sf::st_coordinates(sf::st_as_sf(
        as.data.frame(spatialdataR::point(sd, "coordinates"))))
    expect_equal(unname(xy[, "X"]), c(5, 6, 7, 8))
    expect_equal(unname(xy[, "Y"]), c(1, 2, 3, 4))

    bad <- .mase_with(points = PointsLayerList(
        coordinates = DataFrame(label = letters[1:4], instance_id = letters[1:4])))
    expect_error(as(bad, "SpatialData"), "coordinates.*needs numeric coordinate columns")
})

test_that("a spatialMap that covers only some columns keeps those instance ids", {
    skip_if_not_installed("spatialdataR")
    expect_true(.has_spatialdata_coercion(),
                info = "SpatialData coercion not registered")
    sm <- DataFrame(assay = factor(rep("t1", 3), "t1"), colname = c("a", "b", "d"),
                    element_type = "points", region = "p",
                    instance_id = c("i1", "i2", "i4"))
    mase <- .mase_with(
        points = PointsLayerList(p = DataFrame(x = 1:3, y = 1:3,
                                               instance_id = c("i1", "i2", "i4"))),
        spatialMap = sm)
    expect_warning(sd <- as(mase, "SpatialData"), "1 of 4 columns are not in the spatialMap")
    cd <- SingleCellExperiment::int_colData(spatialdataR::tables(sd)[["t1"]])
    expect_equal(cd$instance_id, c("i1", "i2", NA, "i4"))
    expect_equal(cd$rk, c("p", "p", NA, "p"))

    back <- as(sd, "MultiAssaySpatialExperiment")
    got <- spatialMap(back)
    expect_equal(nrow(got), 3L)
    expect_equal(as.character(got$colname), c("a", "b", "d"))
    expect_equal(as.character(got$instance_id), c("i1", "i2", "i4"))
})

test_that("an element name used in two layers is refused with a clear message", {
    skip_if_not_installed("spatialdataR")
    expect_true(.has_spatialdata_coercion(),
                info = "SpatialData coercion not registered")
    mase <- .mase_with(
        points = PointsLayerList(p = DataFrame(x = 1:2, y = 1:2, instance_id = c("a", "b"))),
        shapes = ShapesLayerList(p = .shapes_df(2)))
    expect_error(as(mase, "SpatialData"), "unique across points, shapes and tables.*p")
})

test_that("the layer replacement methods accept what the getters return", {
    skip_if_not_installed("spatialdataR")
    expect_true(.has_spatialdata_coercion(),
                info = "SpatialData coercion not registered")
    library(spatialdataR)
    mase <- .mase_with(
        points = PointsLayerList(p = DataFrame(x = 1:2, y = 1:2, instance_id = c("a", "b"))),
        shapes = ShapesLayerList(s = .shapes_df(2)),
        images = RasterLayerList(im = matrix(1:4, 2)))
    for (f in c("images", "labels", "points", "shapes", "tables")) {
        x <- mase
        x <- get(paste0(f, "<-"), asNamespace("spatialdataR"))(
            x, value = get(f, asNamespace("spatialdataR"))(x))
        expect_equal(names(get(f, asNamespace("spatialdataR"))(x)),
                     names(get(f, asNamespace("spatialdataR"))(mase)), info = f)
    }
})

test_that("coordinate transformations are reported, but the identity is not", {
    skip_if_not_installed("spatialdataR")
    expect_true(.has_spatialdata_coercion(),
                info = "SpatialData coercion not registered")
    tb <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1:6, 2, 3)),
        colData = DataFrame(row.names = c("a", "b", "c")))
    p <- spatialdataR::SpatialDataPoint(
        data = data.frame(x = 1:3, y = 1:3, instance_id = c("a", "b", "c")))
    moved <- spatialdataR::addCT(p, "scaled", "scale", c(2, 2))
    expect_warning(
        as(spatialdataR::SpatialData(points = list(pp = moved), tables = list(t = tb)),
           "MultiAssaySpatialExperiment"),
        "coordinate transformations are not converted.*pp")
    expect_no_warning(
        as(spatialdataR::SpatialData(points = list(pp = p), tables = list(t = tb)),
           "MultiAssaySpatialExperiment"))
})

test_that("a table keeps its own region and instance keys across a round trip", {
    skip_if_not_installed("spatialdataR")
    expect_true(.has_spatialdata_coercion(),
                info = "SpatialData coercion not registered")
    library(SingleCellExperiment)
    tb <- SingleCellExperiment(list(counts = matrix(1:6, 2, 3)),
        colData = DataFrame(myregion = rep("pp", 3), iid = c("a", "b", "c"),
                            row.names = c("a", "b", "c")))
    int_metadata(tb)$spatialdata_attrs <- list(
        region = "pp", region_key = "myregion", instance_key = "iid")
    p <- spatialdataR::SpatialDataPoint(
        data = data.frame(x = 1:3, y = 1:3, instance_id = c("a", "b", "c")))
    sd <- spatialdataR::SpatialData(points = list(pp = p), tables = list(t = tb))
    back <- as(as(sd, "MultiAssaySpatialExperiment"), "SpatialData")
    out <- spatialdataR::tables(back)[["t"]]
    attrs <- int_metadata(out)$spatialdata_attrs
    expect_equal(attrs$region_key, "myregion")
    expect_equal(attrs$instance_key, "iid")
    expect_false("instance_id" %in% colnames(SummarizedExperiment::colData(out)))
})

test_that("a points element links through the table's instance_key, whatever its name", {
    skip_if_not_installed("spatialdataR")
    library(SingleCellExperiment)
    sce <- SingleCellExperiment(
        list(counts = matrix(1:12, 3, 4)),
        colData = DataFrame(rk = rep("cells", 4),
                            cell_id = paste0("c", 1:4),
                            row.names = paste0("cell", 1:4)))
    int_metadata(sce)[["spatialdata_attrs"]] <- list(
        region = "cells", region_key = "rk", instance_key = "cell_id")
    sd <- spatialdataR::SpatialData(
        images = list(), labels = list(), shapes = list(),
        points = list(cells = spatialdataR::SpatialDataPoint(data = data.frame(
            x = 1:4, y = 2:5, cell_id = paste0("c", 1:4)))),
        tables = list(assay1 = sce))

    mase <- as(sd, "MultiAssaySpatialExperiment")
    pts <- spatialPoints(mase)[["cells"]]
    expect_equal(as.character(pts$instance_id), paste0("c", 1:4))
    smap <- spatialMap(mase)
    expect_equal(as.character(smap$instance_id), paste0("c", 1:4))
    expect_true(all(smap$instance_id %in% pts$instance_id))
})

.poly_layer <- function(crs = NA) {
    g <- sf::st_sfc(lapply(1:4, function(i) sf::st_polygon(list(
        rbind(c(i, 0), c(i + 1, 0), c(i + 1, 1), c(i, 0))))))
    if (!is.na(crs)) sf::st_crs(g) <- crs
    d <- DataFrame(instance_id = paste0("cell", 1:4), val = c(1.5, 2, 3, 4))
    d$geometry <- g
    d
}

test_that("a shapes layer linked to a table keeps all its rows and its CRS", {
    skip_if_not_installed("spatialdataR")
    shp <- .poly_layer(4326)
    sm <- DataFrame(assay = factor(rep("t1", 4), "t1"),
                    colname = c("a", "b", "c", "d"), element_type = "shapes",
                    region = "reg", instance_id = paste0("cell", 1:4))
    mase <- .mase_with(shapes = ShapesLayerList(reg = shp), spatialMap = sm)

    sd <- as(mase, "SpatialData")
    el <- spatialdataR::shape(sd, "reg")
    expect_identical(length(el), 4L)
    back <- as(sd, "MultiAssaySpatialExperiment")
    got <- spatialShapes(back)[["reg"]]
    expect_identical(nrow(got), 4L)
    expect_setequal(as.character(got$instance_id), paste0("cell", 1:4))
    expect_equal(sf::st_crs(got$geometry)$epsg, 4326L)
    expect_equal(nrow(spatialMap(back)), 4L)
})

test_that("a shape the tables do not list is reported, not dropped silently", {
    skip_if_not_installed("spatialdataR")
    shp <- .poly_layer()
    sm <- DataFrame(assay = factor(rep("t1", 3), "t1"),
                    colname = c("a", "b", "c"), element_type = "shapes",
                    region = "reg", instance_id = paste0("cell", 1:3))
    mase <- .mase_with(shapes = ShapesLayerList(reg = shp), spatialMap = sm)
    expect_warning(as(mase, "SpatialData"), "3 of its 4 rows|has 3 of its 4")
})

test_that("a conversion that fails leaves no tables behind", {
    skip_if_not_installed("spatialdataR")
    conn <- MultiAssaySpatialExperiment:::.sd_conn()
    before <- DBI::dbListTables(conn)
    mase <- .mase_with(
        points = PointsLayerList(p = DataFrame(x = 1:2, y = 1:2,
                                               instance_id = c("a", "b"))),
        shapes = ShapesLayerList(bad = DataFrame(a = 1:3)))
    expect_error(as(mase, "SpatialData"), "no geometry column")
    expect_setequal(DBI::dbListTables(conn), before)
})

test_that("element lookups fail the same way for every accessor", {
    skip_if_not_installed("spatialdataR")
    mase <- .mase_with(
        points = PointsLayerList(p = DataFrame(x = 1:2, y = 1:2,
                                               instance_id = c("a", "b"))))
    expect_error(spatialdataR::point(mase, "zz"), "invalid 'i'")
    expect_error(spatialdataR::point(mase, 5), "invalid 'i'")
    expect_error(spatialdataR::table(mase, "zz"), "invalid 'i'")
    expect_null(spatialdataR::shape(mase, 1))
    expect_s4_class(spatialdataR::point(mase, "p"), "DFrame")
})

test_that("replacement methods refuse what they cannot place", {
    skip_if_not_installed("spatialdataR")
    mase <- .mase_with(
        points = PointsLayerList(p = DataFrame(x = 1:2, y = 1:2,
                                               instance_id = c("a", "b"))))
    v <- DataFrame(x = 3:4, y = 3:4, instance_id = c("c", "d"))
    expect_error(spatialdataR::point(mase, 10) <- v, "add a new one by name")
    spatialdataR::point(mase, "q") <- v
    expect_identical(spatialdataR::pointNames(mase), c("p", "q"))
    spatialdataR::element(mase, 1) <- v
    expect_equal(as.character(spatialdataR::point(mase, "p")$instance_id),
                 c("c", "d"))
    spatialdataR::points(mase) <- NULL
    expect_length(spatialdataR::points(mase), 0L)
    spatialdataR::images(mase) <- NULL
    expect_length(spatialdataR::images(mase), 0L)
})

test_that("renaming an element renames the region the spatialMap names", {
    skip_if_not_installed("spatialdataR")
    sm <- DataFrame(assay = factor(rep("t1", 2), "t1"), colname = c("a", "b"),
                    element_type = "points", region = "p",
                    instance_id = c("i1", "i2"))
    mase <- .mase_with(
        points = PointsLayerList(p = DataFrame(x = 1:2, y = 1:2,
                                               instance_id = c("i1", "i2"))),
        spatialMap = sm)
    spatialdataR::pointNames(mase) <- "cells"
    expect_identical(spatialdataR::pointNames(mase), "cells")
    expect_identical(unique(as.character(spatialMap(mase)$region)), "cells")
    expect_error(spatialdataR::pointNames(mase) <- "t1", "unique across")
    expect_error(spatialdataR::tableNames(mase) <- "x", "cannot be renamed")
})
