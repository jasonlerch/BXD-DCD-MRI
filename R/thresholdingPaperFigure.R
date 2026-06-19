dbabl6comp <- function(x, brainbehaviourNoOut) {
  tmp <- brainbehaviourNoOut %>% 
    mutate(v = x$normVolumes) %>%
    filter(Genotype %in% c("C57BL_6J", "DBA_2J")) %>%
    mutate(v = scale(v))
  l <- lm(v ~ Genotype, tmp)
  tidy(l, conf.int=T) %>% 
    filter(term != "(Intercept)") %>%
    select(-term)
  
}

# not related to the BXD paper, but instead for Paul Taylor's thresholding paper
thresholdingPaperFigure <- function(brainbehaviourNoOut, hvolsNoOut, DSURQEanatVol, DSURQElabelVol) {
  ## let's compute some Bl6 vs DBA differences
  bl6vsdba <- Clone(hvolsNoOut)
  bl6vsdba$Do(function(x) {
    tmp <- dbabl6comp(x, brainbehaviourNoOut)
    for (n in colnames(tmp)) {
      x[[n]] <- tmp[n]
    }
  })
  
  sliceList <- map(seq(100, 350, by=50), ~ c(.x, 2))
  sliceList[[length(sliceList)+1]] <- c(145, 1)
  
  slices <- MRIcrotome(DSURQEanatVol * (DSURQElabelVol > 0.5), DSURQElabelVol, bl6vsdba, 
                       sliceList, sliceOffset = 5, assembleDir = "X")
  
  roiCols <- ToDataFrameTable(bl6vsdba, "name", "conf.low", "estimate", "conf.high", "p.value", "statistic")
  high <- round(max(abs(roiCols$estimate)), digits = 2)
  low <- round(min(abs(roiCols$estimate[roiCols$p.value<0.05])), digits = 2)
  
  roi <- inner_join(slices$data, roiCols, by=c("region"="name"))
  roiL <- pivot_longer(roi, cols=conf.low:conf.high) %>% 
    mutate(name=fct_relevel(name, "conf.low", "estimate", "conf.high"))
  # MRIcroscope(slices) %>%
  #   add_anatomy() + 
  #   #geom_sf(data=roiL, colour=NA, aes(fill=value, alpha = factor(p.value<0.05))) + 
  #   geom_sf(data=roiL, colour=NA, aes(fill=value)) + 
  #   #scale_fill_gradient2("Effect size", low=muted("blue"), high=muted("red")) + 
  #   scale_fill_posneg2("Effect size", high=high, low=low, lowalpha = 0.2) + 
  #   #scale_alpha_manual("Significant?", values=c(0.7, 1)) + 
  #   facet_wrap(~name, ncol=1, strip.position = "left") + 
  #   geom_sf(data=roiL, fill=NA, aes(colour=factor(p.value<0.05))) + 
  #   scale_color_manual("Significant?", values=c(NA, "black"), na.value="transparent", guide="none")
  # MRIcroscope(slices) %>%
  #   add_anatomy() + 
  #   geom_sf(data=roiL, colour=NA, aes(fill=value,alpha=abs(value/0.72)^2)) + 
  #   scale_fill_gradientn("Effect size", colours =c("turquoise1", 
  #                                                  "blue", 
  #                                                  "transparent", 
  #                                                  "red", 
  #                                                  "yellow"), breaks=c(-2, -0.72, 0, 0.72, 2)) + 
  #   scale_alpha(guide="none") + 
  #   geom_sf(data=roi, fill=NA, aes(colour=factor(p.value<0.05))) + 
  #   facet_wrap(~name, ncol=1, strip.position="left") + 
  #   scale_color_manual("Thresholding", values=c(NA, "black"), labels=c("q>0.05", "q<0.05"), na.value="transparent")
  # MRIcroscope(slices) %>%
  #   add_anatomy() + 
  #   geom_sf(data=roi, colour=NA, aes(fill=statistic,alpha=I(abs(statistic/2)^2))) + 
  #   scale_fill_gradientn("t-statistic", colours =c("turquoise1", 
  #                                                  "blue", 
  #                                                  "transparent", 
  #                                                  "red", 
  #                                                  "yellow"), 
  #                        breaks=c(-5, -2, 0, 2, 5), 
  #                        values = rescale(c(-5,-0.001, 0, 0.001, 5)),
  #                        limits=c(-5,5), 
  #                        na.value="transparent",
  #                        oob=squish) + 
  #   scale_alpha("t-statistic", range=c(0.5,1)) + 
  #   geom_sf(data=roi, fill=NA, aes(colour=factor(p.value<0.05))) + 
  #   #facet_wrap(~name, ncol=1, strip.position="left") + 
  #   scale_color_manual("Thresholding", values=c(alpha("black", 0.1), "black"), labels=c("q>0.05", "q<0.05"), na.value="transparent")
  MRIcroscope(slices) +
    #add_anatomy() + 
    geom_sf_interactive(data=roi, linewidth=0.25,
                        aes(fill=statistic, colour=factor(p.value<0.05),
                            tooltip=paste0(region, ": ", round(statistic, 2)))) + #,alpha=I(abs(statistic/2)^2))) + 
    scale_fill_posneg2("t-statistic", low=2, high=5) + 
    #scale_alpha("t-statistic", range=c(0.5,1)) + 
    #geom_sf(data=roi, fill=NA, aes(colour=factor(p.value<0.05))) + 
    #facet_wrap(~name, ncol=1, strip.position="left") + 
    scale_color_manual("Thresholding", values=c(alpha("black", 0.1), "black"), labels=c("q>0.05", "q<0.05"), na.value="transparent") + 
    theme_void()
}
