### =========================================================================
### Package load / unload
### -------------------------------------------------------------------------
###
### .onLoad registers SFE methods/coercions when those
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
}
