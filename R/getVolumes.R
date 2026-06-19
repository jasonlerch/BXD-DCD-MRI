# saved rds file has volumes for all 202 brains
# but here we will exclude those without a known Strain
getVolumes <- function(rdsfile, fullbrainbehaviour) {
  hvols <- readRDS(rdsfile)
  hvols$Do(function(x) {
    x$volumes <- x$volumes[fullbrainbehaviour$Strain != ""]
  })
  return(hvols)
}