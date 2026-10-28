relocate_bpcells <- function(x, root) {
  root <- normalizePath(root, mustWork=TRUE)
  changes <- list()
  walk <- function(z) {
    if (isS4(z)) {
      if (inherits(z, 'MatrixDir')) {
        old <- slot(z, 'dir')
        stopifnot(length(old)==1L)
        target <- file.path(root, basename(old))
        if (!dir.exists(target)) stop('Missing BPCells directory: ', target)
        probe <- BPCells::open_matrix_dir(target)
        if (!identical(dim(probe), dim(z))) stop('Dimension mismatch: ', target)
        slot(z, 'dir') <- normalizePath(target, mustWork=TRUE)
        changes[[length(changes)+1L]] <<- data.frame(old=old,new=slot(z,'dir'))
      } else {
        for (n in slotNames(z)) {
          v <- slot(z,n)
          if (isS4(v)||is.list(v)) slot(z,n) <- walk(v)
        }
      }
    } else if (is.list(z)) z <- lapply(z,walk)
    z
  }
  result <- walk(x)
  if (!length(changes)) stop('No BPCells MatrixDir objects found')
  list(counts=result, mapping=unique(do.call(rbind,changes)))
}
