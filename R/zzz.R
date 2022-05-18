# Packgage library 
.onLoad <- function(libname, pkgname) {
  library.dynam("survAUC", pkgname, libname);
}

.onUnload <- function (libpath) {
  library.dynam.unload("survAUC", libpath)
} 