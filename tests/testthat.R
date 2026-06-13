## Load SpatialFeatureExperiment before MASE when running via R CMD check so
## SFE coercions/methods register during .onLoad (SFE is Suggested, not Imported).
if (requireNamespace("SpatialFeatureExperiment", quietly = TRUE))
    loadNamespace("SpatialFeatureExperiment")

library(testthat)
library(MultiAssaySpatialExperiment)

test_check("MultiAssaySpatialExperiment")
