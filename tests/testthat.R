## Load optional packages first so MASE registers coercion/methods in .onLoad
if (requireNamespace("SpatialFeatureExperiment", quietly = TRUE))
    library(SpatialFeatureExperiment)

library(testthat)
library(MultiAssaySpatialExperiment)

test_check("MultiAssaySpatialExperiment")
