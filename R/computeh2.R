h2 <- function(f, df, position=1) {
  l <- broom::tidy(anova(lm(f, data=df)))
  l$sumsq[position] / sum(l$sumsq)
}

h22 <- function(x, Strain) {
  l <- broom::tidy(anova(lm(x~Strain)))
  l$sumsq[1] / sum(l$sumsq)
}

# bootstrapped h2 calculation, stratified sampling by strain
h2boot <- function(f, df, position=1, samples=100) {
  map_dbl(1:samples, ~ h2(f, df %>% group_by(Strain) %>% sample_frac(size=1, replace=T) %>% ungroup, position=position))
}

# to test motor vs nomotor heritability
h2permute <- function(f, df, position=1, samples=100) {
  tmp <- map_dbl(1:samples, function(x) {
    tmpdf <- df %>% group_by(Strain) %>% sample_frac(size=1) %>% ungroup
    df$pmotor <- tmpdf$motor
    h2(f, df %>% filter(pmotor=="motor")) - h2(f, df %>% filter(pmotor=="no"))
  })
  mean(abs(h2(f, df %>% filter(motor=="motor")) - h2(f, df %>% filter(motor=="no"))) < abs(tmp))
  
}

computeh2 <- function(hvols, brainbehaviour, samples=100) {
  out <- Clone(hvols)
  
  nodeList <- Traverse(out)
  
  plan(multisession, workers=10)
  mappedList <- future_map(nodeList, function(x) { #out$Do(function(x){
    brainbehaviour$v <- x$volumes
    brainbehaviour$nv <- x$normVolumes
    
    x$absh2 <- h2(v ~ Strain, brainbehaviour)
    x$normh2 <- h2(nv ~ Strain, brainbehaviour)
    
    x$absh2b <- h2boot(v ~ Strain, brainbehaviour)
    x$normh2b <- h2boot(nv ~ Strain, brainbehaviour)
    
    x$absh2motor <- h2(v ~ Strain, brainbehaviour %>% filter(motor == "motor"))
    x$normh2motor <- h2(nv ~ Strain, brainbehaviour %>% filter(motor == "motor"))
    
    x$absh2nomotor <- h2(v ~ Strain, brainbehaviour %>% filter(motor == "no"))
    x$normh2nomotor <- h2(nv ~ Strain, brainbehaviour %>% filter(motor == "no"))
    
    x$absh2motorb <- h2boot(v ~ Strain, brainbehaviour %>% filter(motor == "motor"), samples=samples)
    x$normh2motorb <- h2boot(nv ~ Strain, brainbehaviour %>% filter(motor == "motor"), samples=samples)
    
    x$absh2nomotorb <- h2boot(v ~ Strain, brainbehaviour %>% filter(motor == "no"), samples=samples)
    x$normh2nomotorb <- h2boot(nv ~ Strain, brainbehaviour %>% filter(motor == "no"), samples=samples)
    
    x$h2motordiff <- x$normh2motor - x$normh2nomotor
    x$h2motordiffp <- h2permute(nv ~ Strain, brainbehaviour, samples=samples)
    
    return(x)
    
  }, .options = furrr_options(seed = TRUE))
  
  # in future_map the tree does not appear to be modified, so do that explicitly
  # variables that have to be copied
  copyVars <- names(mappedList[[1]])[ names(mappedList[[1]]) %in% names(nodeList[[1]]) == FALSE]
  purrr::map(1:length(nodeList), function(x){
    purrr::map(copyVars, function(n) {
      nodeList[[x]][[n]] <- mappedList[[x]][[n]]
    })
  })
  
  return(out)
}