### =========================================================================
### Package load / unload
### -------------------------------------------------------------------------
###
### .onLoad registers SFE and SpatialData methods/coercions when those
### packages are already loaded (isNamespaceLoaded).
###
### -------------------------------------------------------------------------

.onLoad <- function(libname, pkgname) {
    ## Register data reader technologies
    .onLoad_register_technologies()
    
    ## Register coercion methods for optional packages
    if (isNamespaceLoaded("SpatialFeatureExperiment")) {
        .register_SFE_coercion()
        .register_SFE_methods()
    }
    if (isNamespaceLoaded("SpatialData")) {
        .register_SpatialData_coercion()
        .register_SpatialData_methods()
    }
}
