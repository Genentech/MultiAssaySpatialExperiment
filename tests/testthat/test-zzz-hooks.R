test_that("optional-package hooks are replaced, not accumulated", {
    pkg <- "maseFakePkgForHooks"
    event <- packageEvent(pkg, "onLoad")
    withr_clean <- function() setHook(event, NULL, action = "replace")
    on.exit(withr_clean())
    user <- function(...) NULL
    setHook(event, user, action = "append")
    for (i in 1:3)
        MultiAssaySpatialExperiment:::.optional_package_hook(pkg, function() NULL)
    hooks <- getHook(event)
    expect_length(hooks, 2L)
    expect_identical(sum(vapply(hooks,
        MultiAssaySpatialExperiment:::.is_mase_hook, logical(1))), 1L)
    MultiAssaySpatialExperiment:::.remove_optional_hooks(pkg)
    hooks <- getHook(event)
    expect_length(hooks, 1L)
    expect_identical(hooks[[1L]], user)
})

test_that("a failing registration warns and names the integration", {
    expect_warning(
        MultiAssaySpatialExperiment:::.register_optional(
            "somePkg", function() stop("boom")),
        "somePkg integration could not be fully registered.*boom")
    expect_no_warning(
        MultiAssaySpatialExperiment:::.register_optional(
            "somePkg", function() NULL))
})

test_that("the installed hook reports failures as a warning", {
    pkg <- "maseFakePkgForHooks2"
    event <- packageEvent(pkg, "onLoad")
    on.exit(setHook(event, NULL, action = "replace"))
    MultiAssaySpatialExperiment:::.optional_package_hook(
        pkg, function() stop("kaput"))
    expect_warning(getHook(event)[[1L]](pkg, ""), "kaput")
})

test_that("reloading an optional package registers the integration again", {
    skip_if_not_installed("SpatialFeatureExperiment")
    skip_if_not_installed("spatialdataR")
    rscript <- file.path(R.home("bin"), "Rscript")
    code <- '
        suppressMessages({
            library(MultiAssaySpatialExperiment)
            library(SpatialFeatureExperiment)
            loadNamespace("spatialdataR")
        })
        has <- function() c(
            sfe = existsMethod("dimGeometry", "MultiAssaySpatialExperiment"),
            sd = existsMethod("images", "MultiAssaySpatialExperiment"))
        before <- has()
        suppressMessages({
            unloadNamespace("SpatialFeatureExperiment")
            loadNamespace("SpatialFeatureExperiment")
            unloadNamespace("spatialdataR")
            loadNamespace("spatialdataR")
        })
        after <- has()
        cat(if (all(before) && all(after)) "RELOAD_OK" else
            paste("RELOAD_FAILED", toString(before), toString(after)), "\n")'
    out <- suppressWarnings(system2(rscript, c("-e", shQuote(code)),
                                    stdout = TRUE, stderr = TRUE))
    expect_match(paste(out, collapse = "\n"), "RELOAD_OK", fixed = TRUE)
})

test_that("unloading the package removes its load hooks", {
    rscript <- file.path(R.home("bin"), "Rscript")
    code <- '
        suppressMessages(library(MultiAssaySpatialExperiment))
        n <- function(p) length(getHook(packageEvent(p, "onLoad")))
        before <- c(n("SpatialFeatureExperiment"), n("spatialdataR"))
        suppressMessages(unloadNamespace("MultiAssaySpatialExperiment"))
        after <- c(n("SpatialFeatureExperiment"), n("spatialdataR"))
        cat(if (all(before == 1L) && all(after == 0L)) "UNLOAD_OK" else
            paste("UNLOAD_FAILED", toString(before), toString(after)), "\\n")'
    out <- suppressWarnings(system2(rscript, c("-e", shQuote(code)),
                                    stdout = TRUE, stderr = TRUE))
    expect_match(paste(out, collapse = "\n"), "UNLOAD_OK", fixed = TRUE)
})
