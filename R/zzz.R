kt <- NULL
.onLoad <- function(libname, pkgname) {
  py_mod <- system.file("python", package = "npiv")
  kt <<- reticulate::import_from_path('kinematics', py_mod, delay_load = TRUE)
}
