# there are some outliers in the data
# this could use some proper QC, but for the moment let's threshold
# criteria will be SD away from group average for brain volume
computeOutlierness <- function(hvols, brainbehaviour, thresh=2.5) {
  brainbehaviour %>% 
    mutate(bv=hvols$volumes) %>% 
    group_by(Strain) %>% 
    mutate(m=mean(bv, na.rm=T), s=sd(bv, na.rm=T)) %>% 
    ungroup() %>% 
    mutate(outlierness = abs( (bv-m) / s),
           isOutlier = is.na(bv) | outlierness > thresh)
}

threshBrainbehaviour <- function(brainbehaviour) {
  brainbehaviour %>% dplyr::filter(isOutlier != TRUE)
}

threshHvols <- function(hvols, brainbehaviour) {
  out <- Clone(hvols)
  out$Do(function(x){
    x$volumes <- x$volumes[brainbehaviour$isOutlier != TRUE]
    x$normVolumes <- x$normVolumes[brainbehaviour$isOutlier != TRUE]
  })
  return(out)
}