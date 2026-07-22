## Generate the minimal mock spatial output directories shipped in
## inst/extdata/, used by the readMERSCOPEMASE() and readCosMxMASE() examples so
## those examples are self-contained and runnable at R CMD check time, without
## requiring real (large) vendor datasets.
##
## Layouts mirror real vendor output directories:
##   merscope_mock/region_0/cell_metadata.csv, cell_boundaries.csv
##   cosmx_mock/cell_metadata.csv, cell_boundaries.csv, expression_matrix.csv
##
## Regenerate from the package root with:
##   Rscript inst/scripts/make-example-mocks.R

## Square cell-boundary polygons (5 vertices each, closed ring)
square <- function(id, x, y, size = 5)
    data.frame(cell_id = id,
               vertex_x = c(x, x + size, x + size, x, x),
               vertex_y = c(y, y, y + size, y + size, y))

## --- Vizgen MERSCOPE ------------------------------------------------------
ms <- file.path("inst", "extdata", "merscope_mock", "region_0")
dir.create(ms, recursive = TRUE, showWarnings = FALSE)
ms_ids <- c("cell1", "cell2", "cell3")
write.csv(
    data.frame(cell_id = ms_ids,
               center_x = c(10, 15, 20),
               center_y = c(20, 25, 30)),
    file.path(ms, "cell_metadata.csv"), row.names = FALSE)
write.csv(
    do.call(rbind, Map(square, ms_ids, c(5, 15, 25), c(10, 20, 30))),
    file.path(ms, "cell_boundaries.csv"), row.names = FALSE)

## --- NanoString CosMx -----------------------------------------------------
cx <- file.path("inst", "extdata", "cosmx_mock")
dir.create(cx, recursive = TRUE, showWarnings = FALSE)
cx_ids <- c("c1", "c2", "c3")
write.csv(
    data.frame(cell_id = cx_ids,
               x_global_px = c(10, 15, 20),
               y_global_px = c(20, 25, 30)),
    file.path(cx, "cell_metadata.csv"), row.names = FALSE)
write.csv(
    do.call(rbind, Map(square, cx_ids, c(5, 15, 25), c(10, 20, 30))),
    file.path(cx, "cell_boundaries.csv"), row.names = FALSE)
expr <- data.frame(Gene = c("GENE1", "GENE2"))
for (id in cx_ids) expr[[id]] <- c(1L, 2L)
write.csv(expr, file.path(cx, "expression_matrix.csv"), row.names = FALSE)
