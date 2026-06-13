### =========================================================================
### Package load / unload
### -------------------------------------------------------------------------
###
### Optional integrations register during .onLoad when their namespace is
### already loaded. SpatialData integration lives in inst/scripts/ (not in
### Suggests yet) and is sourced on demand via .register_SpatialData_all().
###
### -------------------------------------------------------------------------

.register_SFE_all <- function() {
    ns <- asNamespace("MultiAssaySpatialExperiment")
    if (isTRUE(get0(".sfe_registered", envir = ns, inherits = FALSE)))
        return(invisible(NULL))
    .register_SFE_coercion()
    .register_SFE_methods()
    assign(".sfe_registered", TRUE, envir = ns)
    invisible(NULL)
}

.optional_package_hook <- function(pkg, register_fn) {
    do_register <- function() {
        if (!requireNamespace(pkg, quietly = TRUE))
            return(invisible(NULL))
        tryCatch(
            register_fn(),
            error = function(e) {
                if (grepl("locked", conditionMessage(e), fixed = TRUE))
                    return(invisible(NULL))
                stop(e)
            }
        )
    }
    setHook(packageEvent(pkg, "onLoad"), function(...) do_register(),
            action = "append")
    if (isNamespaceLoaded(pkg))
        do_register()
}

.import_SFE_generic <- function(name, package = "SpatialFeatureExperiment") {
    dest <- asNamespace("MultiAssaySpatialExperiment")
    if (exists(name, envir = dest, inherits = FALSE))
        return(invisible(get(name, envir = dest)))
    gen <- getGeneric(name, package = package)
    if (is.null(gen))
        stop(package, " generic '", name, "' not found")
    assign(name, gen, envir = dest)
    invisible(gen)
}

.import_SFE_generics <- function() {
    generics <- c(
        "dimGeometry", "dimGeometries", "dimGeometryNames",
        "annotGeometry", "annotGeometries", "annotGeometryNames",
        "spatialGraphs", "spatialGraph", "spatialGraphNames",
        "bbox", "unit", "localResults", "localResultNames",
        "splitByCol", "findDebrisCells"
    )
    for (g in generics)
        .import_SFE_generic(g)
    .import_SFE_generic("aggregate", package = "stats")
    invisible(NULL)
}

.register_SpatialData_all <- function() {
    pkgname <- "MultiAssaySpatialExperiment"
    ns <- asNamespace(pkgname)
    if (isTRUE(get0(".spatialdata_registered", envir = ns, inherits = FALSE)))
        return(invisible(NULL))
    script <- system.file("scripts", "MultiAssaySpatialExperiment-SpatialData.R",
                          package = pkgname)
    if (!nzchar(script))
        return(invisible(NULL))
    sys.source(script, envir = ns)
    utils::getFromNamespace(".register_SpatialData_activate", pkgname)()
    assign(".spatialdata_registered", TRUE, envir = ns)
    invisible(NULL)
}

.onLoad <- function(libname, pkgname) {
    ## Register data reader technologies
    .onLoad_register_technologies()

    if (isNamespaceLoaded("SpatialFeatureExperiment"))
        .register_SFE_all()
    else
        .optional_package_hook("SpatialFeatureExperiment", .register_SFE_all)

    .optional_package_hook("SpatialData", .register_SpatialData_all)
}
