plasticityOnROI <- function(x, brainbehaviourNoOut) {
  df <- brainbehaviourNoOut %>%
    mutate(roi = x$volumes,
           group = fct_relevel(group, "no motor training"),
           motor = fct_relevel(motor, "no"))
  l1 <- lm(roi ~ bv + motor, df)
  l2 <- lm(roi ~ bv + group, df)
  rbind(tidy(l1) %>% filter(term == "motormotor"), 
        tidy(l2) %>% filter(term %in% c("groupcomplex wheel", "groupskilled reaching"))) %>% 
    select(term, statistic, p.value) %>% 
    pivot_wider(values_from = c(statistic, p.value), names_from = term)
}

computePlasticityEffects <- function(brainbehaviourNoOut, h2vols) {
    out <- Clone(h2vols)
    out$Do(function(x){
      tmp <- plasticityOnROI(x, brainbehaviourNoOut)
      for (n in names(tmp))
        x[[n]] <- tmp[n]
    })
    return(out)
}

plotPlasticityEffects <- function(plasticity, DSURQEanatVol, DSURQElabelVol) {
  sliceList <- map(seq(100, 350, by=50), ~ c(.x, 2))
  sliceList[[length(sliceList)+1]] <- c(145, 1)
  slices <- MRIcrotome(DSURQEanatVol * (DSURQElabelVol > 0.5), DSURQElabelVol, h2vols, 
                       sliceList, sliceOffset = 5)
  tmp <- ToDataFrameTable(plasticity, "name", "statistic_motormotor", "statistic_groupskilled reaching", "statistic_groupcomplex wheel")
  roi <- inner_join(slices$data, tmp, by=c("region"="name"))
  roiL <- pivot_longer(roi, 
                       statistic_motormotor:`statistic_groupcomplex wheel`, 
                       names_to = "test", 
                       values_to = "statistic") %>%
    mutate(test = 
             case_when(test == "statistic_motormotor" ~ "any motor",
                       test == "statistic_groupskilled reaching" ~ "skilled reaching",
                       test == "statistic_groupcomplex wheel" ~ "complex wheel"),
           test = fct_relevel(test, "any motor"))
  MRIcroscope(slices) %>%
    add_anatomy() +
    geom_sf_interactive(data=roiL, 
                        aes(alpha=I(abs(statistic/2.75)^2), 
                            linewidth=ifelse(abs(statistic)>2.75, "q<0.05", "q>0.05"), 
                            colour=ifelse(abs(statistic)>2.75, "q<0.05", "q>0.05"), 
                            fill=statistic, 
                            tooltip=paste0(region, ": ", round(statistic, 2)))) +
    scale_fill_posneg2(low=2.75, high=6) + 
    scale_colour_manual("Thresholding", values=c("black", alpha("black", 0.2))) + 
    scale_linewidth_manual("Thresholding", values=c(0.4, 0.1)) +
    facet_wrap(~test)
  
}

plotPlasticityRoi <- function(brainbehaviourNoOut, h2vols, roi) {
  df <- brainbehaviourNoOut %>%
    mutate(vol = FindNode(h2vols, roi)$volumes,
           res = residuals(lm(vol ~ Genotype + bv)) + mean(vol)) %>%
    filter(group != "non-specific motor training") %>% 
    mutate(group=fct_relevel(group, "no motor training"))
  
  p1 <- ggplot(df) + 
    aes(x=group, y=res, colour=group) + 
    geom_beeswarm(alpha=0.5) + 
    stat_summary(fun.data=mean_cl_boot, geom="errorbar", linewidth=1.5) + 
    ylab(bquote(Residual ~ Volume ~ (mm^3))) + 
    theme_light() + 
    theme(axis.text.x = element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank()) + 
    scale_color_brewer(palette = "Set1") 
  p2 <- p1 + facet_wrap(~Genotype, nrow=1) + 
    theme(axis.text.x = element_blank(), 
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          strip.text.x = element_text(angle=90))
  
  (p1 + ggtitle(roi)) + p2 + plot_layout(widths=c(3,7), guides = "collect")
}
