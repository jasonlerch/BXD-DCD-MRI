library(targets)
library(tarchetypes)
library(future)
library(readxl)
library(conflicted)

# to be run on single machine with multiple processors
plan(multisession)

# load R scripts containing functions for this pipeline
lapply(list.files("./R", full.names = TRUE), source)
# load missing MRIcrotome stuff (should fix this with a proper install)
source("~/src/MRIcrotome/R/MRIcroscope.R")

options(tidyverse.quiet = TRUE)
# needed packages
tar_option_set(packages = c("tidyverse",
                            "conflicted",
                            #"dotenv",
                            "rmarkdown",
                            "RMINC",
                            "MRIcrotome",
                            "forcats",
                            "grid",
                            "ggplot2",
                            "broom",
                            "stringr",
                            "cluster",
                            "cowplot",
                            "htmlTable",
                            #"rlist",
                            "rstanarm",
                            #"partykit",
                            "readxl",
                            #"ggparty",
                            "patchwork",
                            #"clustree",
                            #"GGally",
                            "gridExtra",
                            "data.tree",
                            #"waffle",
                            "ggsci",
                            "ggh4x",
                            #"tidybayes",
                            "future",
                            "furrr", 
                            "terra",
                            "sf",
                            "smoothr",
                            "tidyterra",
                            "scales",
                            "ggnewscale",
                            "ggbeeswarm",
                            "ggiraph",
                            "mixOmics",
                            "ggrepel",
                            "ggforce",
                            "ggtext",
                            "gt"),
               memory = "transient",       # slower, since data has to be reloaded
               garbage_collection = TRUE)  # but uses less memory, making local parallelism possible

# chose which versions of conflicted functions to prefer
conflict_prefer("legend", "MRIcrotome")
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

# the pipeline itself
list(
  # the input filenames
  tar_target(volumesRData, "fullbxdvols.RDS", format="file"),
  tar_target(newcsv, "files_absolute_BXD.csv", format="file"),
  tar_target(oldcsv, "goldowitz_rel.csv", format="file"),
  tar_target(motorlearning, "Motor Learning_13Mar2022.xlsx", format="file"),
  tar_target(motorcontrol, "Motor Co data_13 March 2022.xlsx", format="file"),
  tar_target(motorheritability, "Motor Learning_Heritability_Mar2022 copy.xlsx", format="file"),
  tar_target(motorcoordination, "Data_Motor Cordination paper.xlsx", format="file"),
  tar_target(DSURQEanatFile, "mincs/DSURQE_40micron_average.mnc", format="file"),
  tar_target(DSURQElabelFile, "mincs/DSURQE_40micron_labels.mnc", format="file"),
  tar_target(DSURQEanatVol, mincArray(mincGetVolume(DSURQEanatFile))),
  tar_target(DSURQElabelVol, mincArray(mincGetVolume(DSURQElabelFile))),
  
  # combine the inputs
  tar_target(fullbrainbehaviour, combineMRIbehaviour2025(newcsv, 
                                                         oldcsv, 
                                                         motorheritability,
                                                         motorcoordination)),
  tar_target(brainbehaviour, fullbrainbehaviour %>% filter(Strain != "")),
  
  # read the volumes
  tar_target(hvols, getVolumes(volumesRData, fullbrainbehaviour)),
  
  # create normalized volumes
  tar_target(hvols2, createNormalizedVols(hvols, brainbehaviour)),
  
  # add some rough and ready QC
  tar_target(brainbehaviour2, computeOutlierness(hvols2, brainbehaviour)),
  tar_target(brainbehaviourNoOut, threshBrainbehaviour(brainbehaviour2)),
  tar_target(hvolsNoOut, threshHvols(hvols2, brainbehaviour2)),
  
  # mouse table
  tar_target(mouseNtable, makeMouseTable(brainbehaviourNoOut)),
  
  # compute heritabilities
  tar_target(h2vols, computeh2(hvolsNoOut, brainbehaviourNoOut, samples=5)),
  tar_target(h2volsSym, makeSymmetric(h2vols)),
  
  # heritability figure
  tar_target(h2fig, h2figure(h2vols, brainbehaviourNoOut, DSURQEanatVol, DSURQElabelVol)),
  
  # compare good to poor learners - or DCD-like to non-DCD-like
  tar_target(goodvbadlearnersfig, makeGoodVBadLearnersFig(brainbehaviourNoOut, h2vols, DSURQEanatVol, DSURQElabelVol)),
  tar_target(goodvbadlearnersOut, ggsave("goodvbadlearners.png", goodvbadlearnersfig, height = 5, width=7, scale=1.5, dpi=300)),
  
  # run massively univariate brain-behaviour correlations
  tar_target(allcormice, runAllCorrelations(brainbehaviourNoOut, h2vols, method="mice")),
  tar_target(allcorstrains, runAllCorrelations(brainbehaviourNoOut, h2vols, method="strains")),
  # now create and save correlation plots
  #tar_target(allcormicetplot, giantCorPlot(allcormice, DSURQEanatVol, DSURQElabelVol, statistic="t") +
  #             ggtitle("Brain-behaviour correlation - t statistics", subtitle = "Correlated across all mice")),
  tar_target(allcormicerplot, giantCorPlot(allcormice, DSURQEanatVol, DSURQElabelVol, statistic="r", nrow = 3) +
               ggtitle("Brain-behaviour correlations", subtitle = "Correlated across all mice")),
  #tar_target(allcorstrainstplot, giantCorPlot(allcorstrains, DSURQEanatVol, DSURQElabelVol, statistic="t") +
  #             ggtitle("Brain-behaviour correlation - t statistics", subtitle = "Correlated across strain means")),
  tar_target(allcorstrainsrplot, giantCorPlot(allcorstrains, DSURQEanatVol, DSURQElabelVol, statistic="r", nrow=3) +
               ggtitle("Brain-behaviour correlations", subtitle = "Correlated across strain means")),
  #tar_target(allcormicetplots, ggsave("allcormicetplots.png", allcormicetplot, width=20, height=10, bg="white")),
  tar_target(allcormicerplots, ggsave("allcormicerplots.png", allcormicerplot, width=7, height=7, bg="white", scale=1.5)),
  #tar_target(allcorstrainstplots, ggsave("allcorstrainstplots.png", allcorstrainstplot, width=20, height=10, bg="white")),
  tar_target(allcorstrainsrplots, ggsave("allcorstrainsrplots.png", allcorstrainsrplot, width=7, height=7, bg="white", scale=1.5)),
  
  ### heatmap of correlations on reduced set of brain regions
  #tar_target()
  
  ### CCA stuff
  # compute the CCA
  tar_target(bbrCCA, computerCCA(h2volsSym, brainbehaviourNoOut, method="shrinkage", ncomp=8)),
  tar_target(CCAxplot, CCAbehavPlot(bbrCCA, ncomp=4)),
  tar_target(CCAbrainplot, CCAbrainSlices(bbrCCA, DSURQEanatVol, DSURQElabelVol, h2volsSym)),
  tar_target(CCAcorp, CCAcorplot(bbrCCA, ncomp = 13)),
  tar_target(variates1, variatesPlot(bbrCCA, comps=1:2, brainbehaviourNoOut, showLegend = T)),
  tar_target(variates2, variatesPlot(bbrCCA, comps=3:4, brainbehaviourNoOut, showLegend = F)),
  tar_target(CCAbp1, CCAbiplot(bbrCCA, comps=1:2, h2volsSym, labelRank=3)),
  tar_target(CCAbp2, CCAbiplot(bbrCCA, comps=3:4, h2volsSym, labelRank=3)),
  tar_target(CCAfig, assembleCCAfig(CCAxplot, CCAbrainplot, CCAcorp, CCAbp1, CCAbp2, variates1, variates2, filename="ccashrinkage.png")),
  tar_target(ccacsvs, writeCCAcsvs(bbrCCA, "CCA-shrinkage-"))
  
)