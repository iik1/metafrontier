.onLoad <- function(libname, pkgname) {
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    register_s3_method("ggplot2", "autoplot", "metafrontier")
    register_s3_method("ggplot2", "autoplot", "malmquist_meta")
    register_s3_method("ggplot2", "autoplot", "boot_tgr")
  }
}

register_s3_method <- function(pkg, generic, class, envir = parent.frame()) {
  fun <- get(paste0(generic, ".", class), envir = envir)
  if (isNamespaceLoaded(pkg)) {
    registerS3method(generic, class, fun, envir = asNamespace(pkg))
  }
  # Also register when the package is loaded later
  setHook(
    packageEvent(pkg, "onLoad"),
    function(...) registerS3method(generic, class, fun,
                                   envir = asNamespace(pkg))
  )
}
