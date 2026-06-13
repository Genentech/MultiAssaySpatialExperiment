## devtools::load_all runs before test helpers; preload SFE and reload so
## SFE coercions/methods register during .onLoad (SFE is Suggested, not Imported).
if (requireNamespace("SpatialFeatureExperiment", quietly = TRUE)) {
    loadNamespace("SpatialFeatureExperiment")
    ns <- asNamespace("MultiAssaySpatialExperiment")
    if (!isTRUE(get0(".sfe_registered", envir = ns, inherits = FALSE)) &&
        requireNamespace("pkgload", quietly = TRUE)) {
        pkgload::load_all(".", quiet = TRUE)
    }
}
