allBehaviourCorrelations <- function(x, brainbehaviourNoOut) {
  brainbehaviourNoOut %>% 
    mutate(vol=x$volumes) %>% 
    dplyr::select(Strain, bv, vol, `Accelerated Rotarod:  Motor Performance Improvements`:Velocity) %>%
    dplyr::select(-matches("Weight"), -matches("Male/Female")) %>% 
    pivot_longer(cols=c(-Strain, -bv, -vol)) %>% 
    group_by(name) %>% 
    summarise(t=tidy(lm(vol ~ bv + value))$statistic[3],
              r=cor(residuals(lm(vol ~ bv)), value, use="complete.obs"),
              rabsvol=cor(vol, value, use="complete.obs"))
}

behaviourMatrix <- function(brainbehaviourNoOut) {
  brainbehaviourNoOut %>% 
    dplyr::select(Genotype, `Accelerated Rotarod:  Motor Performance Improvements`:Velocity) %>%
    dplyr::select(-matches("Weight"), -matches("Male/Female")) %>% 
    pivot_longer(cols=c(-Genotype)) %>% 
    group_by(Genotype, name) %>% 
    summarise(v=mean(value, na.rm=T)) %>%
    pivot_wider(names_from = name, values_from = v)
}

strainMeans <- function(x, brainbehaviourNoOut, bv=NULL) {
  brainbehaviourNoOut %>%
    mutate(vols = x$volumes) %>%
    dplyr::select(Genotype, vols) %>%
    {if (!is.null(bv)) mutate(., vols = residuals(lm(vols ~ bv, .))) else .} %>%
    #mutate(vols = residuals(lm(vols ~ bv, .))) %>%
    group_by(Genotype) %>%
    summarise( "{x$name}" := mean(vols, na.rm=T))
}

brainMatrix <- function(h2vols, brainbehaviourNoOut, bv=NULL, filterFun=NULL) {
  nodes <- Traverse(h2vols, filterFun = filterFun)
  purrr::map(nodes, strainMeans, brainbehaviourNoOut, bv) %>%
    reduce(left_join, by="Genotype")
  
}

allBehaviourCorrelationsMeans <- function(x, brainbehaviourNoOut) {
  brainbehaviourNoOut %>% 
    mutate(vol=x$volumes) %>% 
    dplyr::select(Genotype, bv, vol, `Accelerated Rotarod:  Motor Performance Improvements`:Velocity) %>%
    dplyr::select(-matches("Weight"), -matches("Male/Female")) %>% 
    pivot_longer(cols=c(-Genotype, -bv, -vol)) %>% 
    group_by(Genotype, name) %>% 
    summarise(bv=mean(bv), vol=mean(vol), v=mean(value, na.rm=T)) %>%
    group_by(name) %>%
    summarise(t=tidy(lm(vol ~ bv + v))$statistic[3],
              r=cor(residuals(lm(vol ~ bv)), v),
              rabsvol=cor(vol, v))
}

runAllCorrelations <- function(brainbehaviourNoOut, h2vols, method="mice") {
  out <- Clone(h2vols)
  out$Do(function(x) {
    if (method == "mice")
      cors <- allBehaviourCorrelations(x, brainbehaviourNoOut)
    else if (method == "strains")
      cors <- allBehaviourCorrelationsMeans(x, brainbehaviourNoOut)
    else
      stop("Uknown method - must be either mice or strains")
    for (i in 1:nrow(cors)) {
      x[[paste0(cors$name[i], " - r")]] <- cors$r[i]
      x[[paste0(cors$name[i], " - t")]] <- cors$t[i]
      x[[paste0(cors$name[i], " - rabsvol")]] <- cors$rabsvol[i]
    }
  })
  return(out)
}

resamplingUnivariateCorrelations <- function(hdefs, anatTable, brainbehaviourNoOut) {
  
}

giantCorPlot <- function(cors, DSURQEanatVol, DSURQElabelVol, statistic="t", filterGroup = "none", nrow=1, 
                         df=12, fdr=0.1) {
  sliceList <- purrr::map(seq(100, 350, by=50), ~ c(.x, 2))
  sliceList[[length(sliceList)+1]] <- c(145, 1)
  slices <- MRIcrotome(DSURQEanatVol * (DSURQElabelVol > 0.5), DSURQElabelVol, cors, 
                       sliceList, sliceOffset = 5)
  bnames <- names(cors)[1:75]
  # always keep t stats as we need them for fdr
  tnames <- bnames[str_ends(bnames, "t")]
  bnames <- bnames[str_ends(bnames, statistic)]
  
  tmp <- do.call(ToDataFrameTable, c(cors, "name", bnames, tnames))

  roi <- inner_join(slices$data, tmp, by=c("region"="name"))
  roiL <- pivot_longer(roi, 
                       c(ends_with("r"), ends_with("t")), 
                       names_to = c("test", ".value"), 
                       names_sep = " - ")  %>%
    mutate(
      #test = str_remove(test, " - ."),
      `behaviour group` = assignBehaviourGroup(test),
      `behaviour test` = assignBehaviourTest(test),
      test = str_remove(test, ".+:\\s*"),
      test = str_replace(test, "Motor Performance Improvements", "Performance Improvements"),
      test = str_replace(test, "Rate of successful steps taken", "Rate of successful steps"),
      p = pt2(t, df),
      q = p.adjust(p, "fdr")
      )
    
  if (filterGroup == "Motor learning")
    roiL <- roiL %>% filter(`behaviour group` == "Motor learning")
  else if (filterGroup == "Generalized test of motor function")
    roiL <- roiL %>% filter(`behaviour group` == "Generalized test of motor function")
  
  sL <- MRIcroscope(slices) %>%
    add_anatomy() 
  
  if (statistic == "t") {
    sL <- sL +
      geom_sf_interactive(data=roiL %>%
                            mutate(test = str_remove(test, " - t")), 
                          aes(alpha=I(abs(statistic/2.1)^2), 
                              linewidth=ifelse(abs(statistic)>2.1, "p<0.05", "p>0.05"), 
                              colour=ifelse(abs(statistic)>2.1, "p<0.05", "p>0.05"), 
                              fill=statistic, 
                              tooltip=paste0(region, ": ", round(statistic, 2)))) +
      scale_fill_posneg_alpha("t-statistic", low=2.1, high=5) + 
      scale_colour_manual("Thresholding", values=c("black", alpha("black", 0.2))) + 
      scale_linewidth_manual("Thresholding", values=c(0.4, 0.1)) +
      facet_nested_wrap(vars(`behaviour group`, `behaviour test`, test), nest_line = element_line(linetype = 1), nrow=nrow)
  } else if (statistic == "r") {
    sL <- sL + 
      geom_sf_interactive(data=roiL, #%>%
                            #mutate(test = str_remove(test, " - r")),
                          aes(alpha=I(abs(r/0.3)^2),
                          fill=r,
                          linewidth=ifelse(q <= fdr, paste0("q<", fdr), paste0("q>", fdr)),
                          colour=ifelse(q <= fdr, paste0("q<", fdr), paste0("q>", fdr)),
                          tooltip=paste0(region, ": ", round(r, 2)))) +
      scale_fill_posneg_alpha("r", low=0.3, high=1) + 
      scale_linewidth_manual("Thresholding", values=c(0.4, 0.1)) + 
      scale_colour_manual("Thresholding", values=c("black", alpha("black", 0.2))) + 
      #facet_wrap(~test, nrow=2)
      #facet_nested(~ `behaviour test` + test, nest_line = element_line(linetype = 1))
      facet_nested_wrap(vars(`behaviour group`, `behaviour test`, test), nest_line = element_line(linetype = 1), nrow=nrow)
  } else if (statistic == "rabsvol") {
    sL <- sL + 
      geom_sf_interactive(data=roiL %>%
                            mutate(test = str_remove(test, " - rabsvol")),
                          aes(alpha=I(abs(statistic/0.3)^2),
                              fill=statistic,
                              tooltip=paste0(region, ": ", round(statistic, 2)))) +
      scale_fill_posneg_alpha("r (abs. vol)", low=0.3, high=1) + 
      facet_nested_wrap(vars(`behaviour test`, test), nest_line = element_line(linetype = 1), nrow=nrow)
  }
  return(sL + theme_void(9))
}

# Hi Jason:
# Here is the list of motor tasks and phenotypes in the order we have in the first two papers:
#   
#   Generalized test of motor function
# Gait analysis
# Body speed
# Leg combination
# Stance duration
# Swing duration
# Posterior extreme position
# Step cycle
# Duty factor
# Rotarod
# Baseline performance; improvement
# Open field
# Total distanced traveled
# Time spent in the centre
# Time spent in the periphery
# Time spent moving
# Time spent not moving
# Velocity
# 
# Motor learning
# Accelerating rotarod
# Motor performance improvements
# Learning
# Horizontal ladder
# Rate of missteps
# Rate of successful steps taken
# Rate of correction
# Complex wheel
# Latency to fall
# Learning
# Skilled reaching task
# First time success
# Total success
# Learning rate
# 
# These are the primary brain areas of interest:
#   Primary motor cortex -> Y
# Secondary motor cortex -> Y
# Precental gyrus -> N
# Primary sensory cortex -> Y
# Secondary sensory cortex -> Y
# Dorsolateral prefrontal cortex (not sure of the equivalent in the mouse) -> "Frontal cortex: area 3"
# Superior frontal, middle frontal, frontal pole -> "Frontal pole, cerebral cortex"
# Posterior parietal cortex -> "Posterior parietal association areas"
# Cingulum/anterior/middle/posterior cingulate -> "Anterior cingulate area"
# Corticospinal tract -> Y
# Cerebral peduncle -> N
# Anterior and posterior thalamic radiation (is this the cervicothalamic tract?) -> Y
# Thalamus -> Y
# Caudate -> Y
# Internal capsule -> Y
# External capsule -> N
# Superior longitudinal fasciculus (not sure of the equivalent in the mouse) -> N
# Inferior fronto-occipital fasciculus (not sure of the equivalent in the mouse) -> N
# Brainstem (don’t have specific regions labelled in the human) -> Y (subdivide?)
# Cerebellum (all of it, but specifically
#             Crus 1 and 2 -> Y
#             Lobules VII -> Y
#             Lobule VIII -> Y
#             Lobule IX -> Y
#             Cerebellar peduncles – superior, middle
#             Corpus callosum (splenium if you have it)
#             
#             Hopefully this helps to organize the data 😊
#             Thanks!
#               Jill
#
# Globus pallidus, VTA, accumbens, pontine nucleus

selectJillsROIs <- function(allcorstrains) {
  nnn <- names(allcorstrains)
  tmp <- purrr::map_dfc(c("name", "meanVolume", nnn[str_ends(nnn, " - r")]), ~ allcorstrains$Get(.x))
  colnames(tmp) <- c("ROI", "meanVolume", str_remove(nnn[str_ends(nnn, "r")], " - r"))
  longtmp <- tmp %>% filter(
    ROI %in% c("Primary motor area", 
               "Secondary motor area", 
               "Primary somatosensory area", 
               "Supplemental somatosensory area",
               "Frontal cortex: area 3",
               "Frontal pole, cerebral cortex",
               "Posterior parietal association areas",
               "Anterior cingulate area",
               "corticospinal tract",
               "cervicothalamic tract",
               "Thalamus",
               "Caudoputamen",
               "internal capsule",
               "Brain stem",
               "Crus 1",
               "Crus 2",
               "Folium-tuber vermis (VII)",
               "Pyramus (VIII)",
               "Uvula (IX)",
               "Cerebellar nuclei",
               "corpus callosum",
               "superior cerebelar peduncles",
               "middle cerebellar peduncle",
               "inferior cerebellar peduncle",
               "Pallidum, dorsal region",
               "Nucleus accumbens",
               "Midbrain, motor related"
                )) %>% 
    pivot_longer(cols = c(-ROI, -meanVolume), names_to = "behaviour", values_to = "r") %>% 
    mutate(behaviour = fct_reorder(behaviour, r)) %>%
    mutate(`behaviour group` = assignBehaviourGroup(behaviour))
  return(longtmp)
}


# assigns whether a behaviour test is for motor coordination or motor learning
# assumes a long df as input
assignBehaviourGroup <- function(behaviour) {
  return(ifelse(behaviour %in% c("Gait analysis",
                                              "Body speed",
                                              "Leg combination",
                                              "Stance duration",
                                              "Swing duration",
                                              "Posterior extreme position",
                                              "Step cycle",
                                              "Duty factor",
                                              "Rotarod",
                                              "Baseline performance",
                                              "Performance Improvement",
                                              "Open field",
                                              "Distance",
                                              "Center Time",
                                              "Peripheral Time",
                                              "Time Moving",
                                              "Time Not Moving",
                                              "Velocity"), 
                "Generalized test of motor function", "Motor learning"))
}

# assigns the battery a behaviour came from
# assumes a long df as input
assignBehaviourTest <- function(behaviour) {
  mapping <- c(
    "Accelerated Rotarod:  Motor Performance Improvements" = "Accelerated Rotarod",
    "Accelerated Rotarod: Learning" = "Accelerated Rotarod",
    "Horizontal Rung: Rate of successful steps taken" = "Horizontal Rung",
    "Horizontal Rung: Rate of missteps" = "Horizontal Rung",
    "Horizontal Rung: Rate of correction" = "Horizontal Rung",
    "Complex Wheel:  Latency to Fall" = "Complex Wheel",
    "Complex Wheel: Learning" = "Complex Wheel",
    "Skilled Reaching: First Time Success" = "Skilled Reaching",
    "Skilled Reaching: Total Success" = "Skilled Reaching",
    "Skilled Reaching: Learning Rate" = "Skilled Reaching",
    "Body speed" = "Gait",
    "Leg combination" = "Gait",
    "Duty factor" = "Gait",
    "Step cycle" = "Gait",
    "Stance duration" = "Gait",
    "Swing duration" = "Gait",
    "Posterior extreme position" = "Gait",
    "Baseline performance" = "Rotarod",
    "Performance Improvement" = "Rotarod",
    "Distance" = "Open Field",
    "Center Time" = "Open Field",
    "Peripheral Time" = "Open Field",
    "Time Moving" = "Open Field",
    "Time Not Moving" = "Open Field",
    "Velocity" = "Open Field"
  )
  return(mapping[behaviour])
}

plotJillHeatmap <- function(jillData, grouping="group", rows=0, columns=0, addText=F) {
  hm <- jillData %>% 
    heatmap(behaviour, ROI, r, palette_value = colorRamp2::colorRamp2(c(-1,-0.3,0,0.3,1), c("turquoise1", "blue", "white", "red", "yellow")),#c("blue", "white", "red"), 
            row_names_gp = gpar(fontsize = 8),
            column_names_gp = gpar(fontsize = 8)) 
  
  if (rows>0)
    hm <- hm %>% split_rows(rows)
  if (columns>0)
    hm <- hm %>% split_columns(columns)
  
  if (grouping == "group")
    hm <- hm %>% annotation_group(`behaviour group`) 
  else
    hm <- hm %>% annotation_tile(`behaviour group`, palette = RColorBrewer::brewer.pal(3, "Set3")[c(1,3)])
  
  if (addText == TRUE)
    hm <- hm %>% layer_text(.value=round(r,2), .size=6)
    
  hm %>% as_ComplexHeatmap() %>% ComplexHeatmap::draw(heatmap_legend_side = "right")
}

computerCCA <- function(h2volsSym, brainbehaviourNoOut, method="ridge", ncomp=4) {
  # create normalized brain and behaviour matrices
  brainNorm <- brainMatrix(h2volsSym, brainbehaviourNoOut, h2volsSym$volumes, isLeaf)
  brainMM <- brainNorm %>% 
    dplyr::select(-Genotype) %>% 
    as.matrix() %>% 
    scale()
  rownames(brainMM) <- brainNorm$Genotype
  
  behaviourMM <- behaviourMatrix(brainbehaviourNoOut) %>% 
    ungroup() %>% 
    dplyr::select(-Genotype) %>% 
    as.matrix() %>%
    scale()
  rownames(behaviourMM) <- behaviourMatrix(brainbehaviourNoOut)$Genotype
  
  if (method == "ridge") {
    # we'll use ridge (L2) normalization for the CCA. Let's tune it first
    tuneCCA <- tune.rcc(behaviourMM, brainMM, validation="loo")
    bbrCCA <- mixOmics::rcc(behaviourMM, brainMM, ncomp = ncomp, method = "ridge", lambda1=tuneCCA$opt.lambda1, lambda2=tuneCCA$opt.lambda2)
  }
  else if (method == "shrinkage") {
    bbrCCA <- mixOmics::rcc(behaviourMM, brainMM, ncomp = ncomp, method = "shrinkage")
  }
  else {
    stop("Unknown method - must be either ridge or shrinkage")
  }
  return(bbrCCA)
}

CCAbehavPlot <- function(bbrCCA, ncomp=4) {
  finalV <- paste0("V", ncomp)
  cor(bbrCCA$X, bbrCCA$variates$X[,1:ncomp], use="pairwise") %>%
    as.data.frame() %>%
    rownames_to_column(var="behaviour") %>%
    pivot_longer(V1:!!finalV, names_to = "comp", values_to = "r") %>%
    ggplot() + 
    aes(x=r, y=behaviour) + 
    geom_col() + 
    ylab("") + 
    facet_grid(.~comp)
}

CCAbrainSlices <- function(bbrCCA, DSURQEanatVol, DSURQElabelVol, h2volsSym, ncomp=4, rthresh=0.3) {
  sliceList <- purrr::map(seq(100, 350, by=50), ~ c(.x, 2))
  sliceList[[length(sliceList)+1]] <- c(145, 1)
  slices <- MRIcrotome(DSURQEanatVol * (DSURQElabelVol > 0.5), DSURQElabelVol, h2volsSym, 
                       sliceList, sliceOffset = 5)
  
  finalV <- paste0("V", ncomp)
  
  ccaData <- inner_join(slices$data, cor(bbrCCA$Y, bbrCCA$variates$Y[,1:ncomp], use="pairwise") %>% 
                          as.data.frame() %>% 
                          rownames_to_column(), 
                        by=c("region" = "rowname")) %>%
    pivot_longer(V1:!!finalV, names_to = "comp", values_to = "r")
  
  MRIcroscope(slices) %>%
    add_anatomy() +
    geom_sf_interactive(data=ccaData, aes(fill=r, alpha=I(abs(r/rthresh)^2))) + 
    scale_fill_posneg_alpha(low=rthresh, high=1) + 
    facet_grid(.~comp)
}

CCAcorplot <- function(bbrCCA, ncomp=14) {
  bbrCCA$cor %>% 
    as.data.frame() %>% 
    rownames_to_column(var="comp") %>% 
    mutate(comp = as.numeric(comp)) %>% 
    filter(comp <= ncomp) %>%
    ggplot() + 
    aes(x=comp, y=.) + 
    geom_col() + 
    ylab("r") + 
    xlab("component") + 
    scale_x_continuous(limits=c(0, ncomp+1), breaks = seq(1,ncomp, by=2)) + 
    theme_minimal()
}

variatesPlot <- function(bbrCCA, comps=1:2, brainbehaviourNoOut, showLegend=FALSE) {
  df <- brainbehaviourNoOut %>%
    mutate(learners = fct_recode(learners, `non DCD-like`="good", `DCD-like`="poor")) %>%
    group_by(Strain) %>%
    summarise(learners = first(learners))
  
  sCols <- c(
    hcl.colors(12, "Purples")[5],
    hcl.colors(7, "Greens")[3]
  )
  
  p1 <- bbrCCA$variates$X[,comps] %>% 
    as.data.frame() %>% 
    rownames_to_column("Strain") %>% 
    left_join(df) %>%
    ggplot() + 
    aes(x=V1, y=V2, colour=learners) + 
    geom_hline(yintercept = 0) + 
    geom_vline(xintercept = 0) + 
    #geom_point() + 
    geom_label_repel(aes(label=Strain), size=2.5, max.overlaps = 30) + 
    xlab(paste("Comp", comps[1])) + 
    ylab(paste("Comp", comps[2])) + 
    scale_color_manual("", values = sCols) + 
    theme_minimal()
  
  if (showLegend)
    p1 <- p1 + 
    coord_cartesian(clip = "off") + 
    theme(legend.position = "inside", 
          legend.position.inside = c(1, 0.2),
          legend.title = element_blank(),
          legend.box.background = element_rect(colour="black"),
          legend.margin = margin(0.5,0.5,0.5,0.5),
          legend.box.margin = margin(0,0,0,0)) 
  else
    p1 <- p1 + theme(legend.position = "none")
  
  return(p1)
  
}

CCAbiplot <- function(bbrCCA, comps=1:2, h2volsSym, labelRank=3) {
  # bind the data together
  behav <- cor(bbrCCA$X, bbrCCA$variates$X[,comps], use="pairwise") %>%
          as.data.frame() %>%
          rownames_to_column(var="variable") %>%
          mutate(block="behaviour", fontface="plain",
                 color_hex_triplet = "#ffffff",
                 r1p=rank(V1), r1n=rank(-V1),
                 r2p=rank(V2), r2n=rank(-V2)) 
  brain <- cor(bbrCCA$Y, bbrCCA$variates$Y[,comps], use="pairwise") %>%
          as.data.frame() %>%
          rownames_to_column(var="variable") %>%
          mutate(block="brain", fontface="italic",
                 r1p=rank(V1), r1n=rank(-V1),
                 r2p=rank(V2), r2n=rank(-V2)) %>%
    left_join(ToDataFrameTree(h2volsSym, "name", "color_hex_triplet"), by=c("variable" = "name")) %>%
    dplyr::select(-levelName)
  df <- rbind(brain,behav)


  labelDf <- df %>%
    filter(r1p <= labelRank | 
             r2p <= labelRank | 
             r1n <= labelRank | 
             r2n <= labelRank)
  #browser()
  ggplot(df) + 
    aes(x=V1, y=V2, colour=I(color_hex_triplet), shape=block) + 
    geom_point() + 
    xlab(paste("Comp", comps[1])) + 
    ylab(paste("Comp", comps[2])) + 
    geom_text_repel(data=labelDf, aes(x=V1, y=V2, label=variable, fontface=fontface), size=3, inherit.aes = T) + 
    geom_circle(aes(x0=x0, y0=y0,r=r), data=data.frame(x0=0, y0=0, r=c(0.5, 0.8)), colour="white",  inherit.aes = F) +
    geom_hline(yintercept = 0, colour="white", linetype = 2) + 
    geom_vline(xintercept = 0, colour="white", linetype = 2) + 
    coord_fixed(xlim=c(-1,1), ylim=c(-1,1)) + 
    theme_minimal() + 
    theme(panel.background = element_rect(fill="black"), panel.grid = element_blank(), legend.position = "none")
}

assembleCCAfig <- function(CCAxplot, CCAbrainplot, CCAcorp, CCAbp1, CCAbp2, variates1, variates2, filename="ccafig.png") {
  design <- "
  ADDD
  BDDD
  CDDD
  EEEE
  FFGG
  "
  fig <- CCAcorp + variates1 + variates2 + CCAbrainplot + free(CCAxplot) + free(CCAbp1) + free(CCAbp2) + 
    plot_layout(design=design, heights = c(0.35, 0.4, 0.4, 0.8,1.25)) + plot_annotation(tag_levels = "A")
  ggsave(filename, plot = fig, height=8, width=6.5, scale=2)
  return(fig)
}

writeCCAcsvs <- function(bbrCCA, file_stem) {
  bbrCCA$variates$Y %>% as.data.frame.matrix() %>% rownames_to_column(var="strain") %>% write_csv(paste0(file_stem, "yvariates.csv"))
  bbrCCA$variates$X %>% as.data.frame.matrix() %>% rownames_to_column(var="strain") %>% write_csv(paste0(file_stem, "xvariates.csv"))
}

makeSymmetric <- function(h2vols) {
  h2volsSym <- Clone(h2vols)
  Prune(h2volsSym, function(x) !str_starts(x$name, "left") & !str_starts(x$name, "right"))
  return(h2volsSym)

}
