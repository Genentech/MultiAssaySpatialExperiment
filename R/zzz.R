### ============================================================================
### Package load / unload
### ----------------------------------------------------------------------------
###
### Optional integrations register during .onLoad when their namespace is
### already loaded. Both integrations register into their own environment
### (.sfe_ext_env, .sd_ext_env), because the namespace is sealed by then.
###
### ----------------------------------------------------------------------------

.sfe_ext_env <- new.env(parent = asNamespace("MultiAssaySpatialExperiment"))
.sd_ext_env <- new.env(parent = asNamespace("MultiAssaySpatialExperiment"))

## The integration is registered once per loaded namespace of the optional
## package. If that package is unloaded and loaded again its generics are new
## objects, so what was registered on the old ones is gone and the generics
## imported from them are stale: forget both and register afresh.
.needs_registration <- function(env, pkg) {
    ns <- asNamespace(pkg)
    if (identical(get0(".namespace", envir = env, inherits = FALSE), ns))
        return(FALSE)
    rm(list = ls(env, all.names = TRUE), envir = env)
    assign(".namespace", ns, envir = env)
    TRUE
}

.register_SFE_all <- function() {
    if (!.needs_registration(.sfe_ext_env, "SpatialFeatureExperiment"))
        return(invisible(NULL))
    ## unset again if registration fails, so that a later load retries
    ok <- FALSE
    on.exit(if (!ok) rm(".namespace", envir = .sfe_ext_env))
    .register_SFE_coercion()
    .register_SFE_methods()
    ok <- TRUE
    invisible(NULL)
}

.optional_packages <- c("SpatialFeatureExperiment", "spatialdataR")

## Hooks this package installs carry an attribute so that they can be found
## again: reloading the package replaces them instead of stacking duplicates,
## and .onUnload removes them without touching other packages' hooks.
.is_mase_hook <- function(f)
    isTRUE(attr(f, "MultiAssaySpatialExperiment"))

.remove_optional_hooks <- function(pkg) {
    event <- packageEvent(pkg, "onLoad")
    hooks <- getHook(event)
    mine <- vapply(hooks, .is_mase_hook, logical(1))
    if (any(mine)) {
        setHook(event, NULL, action = "replace")
        for (f in hooks[!mine])
            setHook(event, f, action = "append")
    }
    invisible(NULL)
}

## A failed registration must not surface as a bare hook error. It leaves
## the integration off (the registered flag is only set on success, so a
## later load can retry) and says which integration and why.
.register_optional <- function(pkg, register_fn) {
    tryCatch(register_fn(), error = function(e)
        warning("MultiAssaySpatialExperiment: the ", pkg, " integration ",
                "could not be fully registered, so some of its methods may ",
                "be missing: ", conditionMessage(e), call. = FALSE))
    invisible(NULL)
}

.optional_package_hook <- function(pkg, register_fn) {
    .remove_optional_hooks(pkg)
    hook <- function(...) .register_optional(pkg, register_fn)
    attr(hook, "MultiAssaySpatialExperiment") <- TRUE
    setHook(packageEvent(pkg, "onLoad"), hook, action = "append")
    if (isNamespaceLoaded(pkg))
        .register_optional(pkg, register_fn)
    invisible(NULL)
}

.import_SFE_generic <- function(name, package = "SpatialFeatureExperiment") {
    dest <- .sfe_ext_env
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
    if (!.needs_registration(.sd_ext_env, "spatialdataR"))
        return(invisible(NULL))
    ok <- FALSE
    on.exit(if (!ok) rm(".namespace", envir = .sd_ext_env))
    .import_SD_generics()
    .register_SpatialData_coercion()
    .register_SpatialData_methods()
    ok <- TRUE
    invisible(NULL)
}

.onLoad <- function(libname, pkgname) {
    ## Register data reader technologies
    .onLoad_register_technologies()

    .optional_package_hook("SpatialFeatureExperiment", .register_SFE_all)
    .optional_package_hook("spatialdataR", .register_SpatialData_all)
}

.onUnload <- function(libpath) {
    for (pkg in .optional_packages)
        .remove_optional_hooks(pkg)
    invisible(NULL)
}
