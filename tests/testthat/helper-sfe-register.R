## Preload SFE so its packageEvent("onLoad") hook fires and SFE
## coercions/methods register (SFE is Suggested, not Imported).
## Registration targets .sfe_ext_env (see zzz.R), not this namespace,
## so it works regardless of whether this namespace was sealed before
## SFE loaded.
if (requireNamespace("SpatialFeatureExperiment", quietly = TRUE)) {
    loadNamespace("SpatialFeatureExperiment")
}
