createNormalizedVols <- function(hvols, brainbehaviour) {
  out <- Clone(hvols)
  brainbehaviour$bv <- out$volumes
  out$Do(function(x) {
    brainbehaviour$v <- x$volumes
    x$normVolumes <- residuals(lm(v ~ bv, brainbehaviour))
  })
  return(out)
}