makeGoodVBadLearnersFig <- function(brainbehaviourNoOut, h2vols, DSURQEanatVol, DSURQElabelVol) {
  
  theme_set(theme_light())

  # set the levels so that colouring and sorting works as desired
  brainbehaviourNoOut <- brainbehaviourNoOut %>%
    mutate(learners = fct_recode(learners, `non-DCD-like`="good", `DCD-like`="poor"),
           brainVolume=h2vols$volumes) %>%
    group_by(learners) %>%
    mutate(Strain = fct_reorder(Strain, brainVolume)) %>%
    ungroup()
  
  
  sCols <- c(
    hcl.colors(12, "Purples")[1:9],
    hcl.colors(7, "Greens")[1:5]
  )
  
  # plot of overall brain volume
  bvBoxPlot <- brainbehaviourNoOut %>%
    ggplot() + 
    aes(x=learners, y=brainVolume) + 
    geom_boxplot() + 
    ylab(bquote(Brain ~ Volume ~ (mm^3))) 
  
  bvStrainBoxPlot <- brainbehaviourNoOut %>%
    ggplot() + 
    aes(x=Strain, y=brainVolume, colour=Strain) + 
    geom_boxplot() + 
    ylab(bquote(Brain ~ Volume ~ (mm^3))) + 
    facet_grid(.~learners, scales = "free_x", space="free_x") + 
    xlab("Strain") + 
    scale_colour_manual(values=sCols, guide="none") +
    theme(axis.text.x = element_text(angle=45, hjust=1))

  # ROIwise analyses of learners
  hlm <- hanatLm(~ brainVolume + learners, brainbehaviourNoOut, h2vols)
  hlmRaw <- hanatLm(~ learners, brainbehaviourNoOut, h2vols)
  #browser()
  sliceList <- purrr::map(seq(100, 350, by=50), ~ c(.x, 2))
  sliceList[[length(sliceList)+1]] <- c(145, 1)
  slices <- MRIcrotome(DSURQEanatVol * (DSURQElabelVol > 0.5), DSURQElabelVol, h2vols, 
                       sliceList, sliceOffset = 5, assembleDir = "X")
  #browser()
  statData <- ToDataFrameTable(hlm, "name", "tvalue.learnersDCD.like") %>%
    mutate(pvals = pt2(tvalue.learnersDCD.like, 190),
           qvals = p.adjust(pvals, "fdr"), 
           type="normalized") %>%
    right_join(slices$data, by=c("name" = "region"))
  statDataRaw <- ToDataFrameTable(hlmRaw, "name", "tvalue.learnersDCD.like") %>%
    mutate(pvals = pt2(tvalue.learnersDCD.like, 190),
           qvals = p.adjust(pvals, "fdr"),
           type="absolute") %>% 
    right_join(slices$data, by=c("name" = "region"))
  #browser()
  statDataBoth <- rbind(statData, statDataRaw)

  learnersBrainPlot <- MRIcroscope(slices) %>%
    add_anatomy() +
    geom_sf_interactive(data=statDataBoth, aes(geometry=geometry, 
                                           fill=tvalue.learnersDCD.like, 
                                           alpha=I(abs(tvalue.learnersDCD.like/2)^2),
                                           linewidth=ifelse(qvals<=0.05, "q<0.05", "q>0.05"), 
                                           colour=ifelse(qvals<=0.05, "q<0.05", "q>0.05"))) + 
    scale_fill_posneg_alpha("t-statistic", low=2, high=5) + 
    scale_colour_manual("Thresholding", values=c("black", alpha("black", 0.2))) + 
    scale_linewidth_manual("Thresholding", values=c(0.4, 0.1)) +
    facet_grid(type ~ .) + 
    theme_void() + 
    theme(legend.position="bottom")
  
  # some brain behaviour correlations
  brainbehaviourNoOut <- brainbehaviourNoOut %>%
    mutate(PMA = FindNode(h2vols, "Primary motor area")$volumes,
           PMA = residuals(lm(PMA~brainVolume)),
           CC = FindNode(h2vols, "corpus callosum")$volumes,
           CC = residuals(lm(CC ~ brainVolume)),
           thalamus = FindNode(h2vols, "corpus callosum")$volumes,
           thalamus = residuals(lm(thalamus ~ brainVolume)))
  
  bbnoStrain <- brainbehaviourNoOut %>%
    group_by(Strain) %>%
    summarise(across(where(is.double), list(mean=function(x) mean(x, na.rm=TRUE), 
                                            sem=function(x) sd(x, na.rm=TRUE)/sqrt(n()),
                                            lCI=function(x) mean_cl_boot(x)$ymin,
                                            uCI=function(x) mean_cl_boot(x)$ymax)))
  stepPlot <- brainbehaviourNoOut %>%
    ggplot() + 
    aes(x=`Step cycle`, y=PMA, colour=Strain) + 
    scale_colour_manual(values=sCols, guide="none") + 
    geom_point(alpha=0.5) + 
    geom_smooth(method="lm", aes(group=NA), linetype=2, fullrange=T, se=F, colour="black") + 
    ylab("Primary motor area volume<br>(Residualized mm<sup>3</sup>)") + 
    geom_point(data=bbnoStrain, aes(x=`Step cycle_mean`, y=PMA_mean), shape=17, size=3) + 
    geom_errorbar(data=bbnoStrain, aes(x=`Step cycle_mean`, ymin=PMA_lCI, ymax=PMA_uCI, y=PMA_mean)) + 
    geom_errorbar(data=bbnoStrain, aes(x=`Step cycle_mean`, xmin=`Step cycle_lCI`, xmax=`Step cycle_uCI`, y=PMA_mean)) + 
    geom_smooth(data=bbnoStrain, aes(x=`Step cycle_mean`, y=PMA_mean, group=NA), method="lm", se=F, fullrange=T, colour="black") + 
    theme(axis.title = element_markdown())
  
  complexwheelPlot <- brainbehaviourNoOut %>%
    ggplot() + 
    aes(x=`Complex Wheel: Learning`, y=thalamus, colour=Strain) + 
    scale_colour_manual(values=sCols, guide="none") + 
    geom_point(alpha=0.5) + 
    geom_smooth(method="lm", aes(group=NA), linetype=2, fullrange=T, se=F, colour="black") + 
    ylab("Thalamus volume<br>(Residualized mm<sup>3</sup>)") +  
    geom_point(data=bbnoStrain, aes(x=`Complex Wheel: Learning_mean`, y=thalamus_mean), shape=17, size=3) + 
    geom_errorbar(data=bbnoStrain, aes(x=`Complex Wheel: Learning_mean`, ymin=thalamus_lCI, ymax=thalamus_uCI, y=thalamus_mean)) + 
    geom_errorbar(data=bbnoStrain, aes(x=`Complex Wheel: Learning_mean`, xmin=`Complex Wheel: Learning_lCI`, xmax=`Complex Wheel: Learning_uCI`, y=thalamus_mean)) + 
    geom_smooth(data=bbnoStrain, aes(x=`Complex Wheel: Learning_mean`, y=thalamus_mean, group=NA), method="lm", se=F, fullrange=T, colour="black") + 
    theme(axis.title = element_markdown())
    
  design <- 
    "
    ABBB
    CCCC
    DDEE
  "
    
  figure <- bvBoxPlot + bvStrainBoxPlot + learnersBrainPlot + stepPlot + complexwheelPlot + 
    plot_layout(axes = "collect", design = design, heights = c(1,1.5,1)) + 
    plot_annotation(tag_levels = "A")
  
  return(figure)
}