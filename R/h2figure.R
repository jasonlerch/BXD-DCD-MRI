h2figure <- function(h2vols, brainbehaviour, DSURQEanatVol, DSURQElabelVol) {
  sliceList <- purrr::map(seq(100, 350, by=50), ~ c(.x, 2))
  sliceList[[length(sliceList)+1]] <- c(145, 1)
  slices <- MRIcrotome(DSURQEanatVol * (DSURQElabelVol > 0.5), DSURQElabelVol, h2vols, 
                       sliceList, sliceOffset = 5)
  tmp <- ToDataFrameTable(h2vols, "normh2", "name", "absh2", "blah")
  pSlices <- ggplot() + 
    geom_spatraster(data=slices$anatomy) + 
    scale_fill_gradient("anatomy", low="black", high="white", limits=c(600, 1400), na.value = "transparent", guide = "none") + 
    new_scale_fill() + 
    geom_spatvector(data=slices$data, fill=NA) + 
    #scale_fill_manual(values=slices$palette, na.value="transparent", guide="none") + 
    new_scale_fill() + 
    geom_spatvector(data=inner_join(slices$data, tmp, by=c("region" = "name")) %>% 
                      pivot_longer(cols=c(normh2, absh2, blah), values_to = "h2", names_to = "type") %>% 
                      mutate(type=fct_relevel(type, "blah", "absh2", "normh2"), 
                             type=fct_recode(type, Anatomy = "blah", Absolute = "absh2", Normalized="normh2")), 
                    aes(fill=h2)) + scale_fill_continuous(na.value="transparent") + 
    facet_grid(~type) + 
    ggtitle("B) heritability of MRI derived volumes") +
    coord_sf() + 
    theme(axis.ticks = element_blank(), panel.grid = element_blank(), 
          panel.background = element_blank(), axis.text = element_blank())
  
  volData <- brainbehaviour %>% 
    mutate(`Cerebral cortex` = FindNode(h2vols, "Cerebral cortex")$volumes,
           Hindbrain = FindNode(h2vols, "Hindbrain")$volumes,
           Cerebellum = FindNode(h2vols, "Cerebellum")$volumes) %>% 
    pivot_longer(cols=`Cerebral cortex`:Cerebellum, names_to = "Brain region", values_to = "volume")
  
  h2s <- volData %>% group_by(`Brain region`) %>% summarise(h2=h22(volume, Strain))
  
  pBV <- volData %>%
    mutate(`Brain region` = fct_inorder(`Brain region`)) %>%
    ggplot() + 
    aes(x=Genotype, y=volume) + 
    geom_beeswarm(alpha=0.4) + 
    stat_summary(fun.data=mean_cl_boot, geom="errorbar", colour="dark blue") + 
    stat_summary(fun.data=mean_cl_boot, geom="point", size=2, colour="dark blue") + 
    ylab(bquote(Volume ~ (mm^3))) + 
    geom_text(data=h2s, aes(x=Inf, y=-Inf, label=paste0("h2=", round(h2, 3))), hjust=1.1, vjust=-0.5) + 
    facet_grid(`Brain region` ~ ., scales="free_y") + 
    ggtitle("A) brain volumes of BXD panel") + 
    xlab("") + 
    theme_classic() + 
    theme(axis.text.x = element_text(angle=22.5, hjust=1))
  
  # pH2 <- brainbehaviour %>%
  #   #select(Strain, `Skilled Reaching: First Time Success`, `Accelerated Rotarod: Learning`, 
  #   #       `Complex Wheel: Learning`, `Body speed`, `Center Time`, `Negative geotaxis`, Velocity, `Step cycle`) %>%
  #   select(Strain, `Accelerated Rotarod:  Motor Performance Improvements`:`Skilled Reaching: Learning Rate`) %>%
  #   mutate(across(-Strain, as.numeric)) %>%
  #   mutate(`Cerebral cortex` = FindNode(h2vols, "Cerebral cortex")$volumes,
  #          Hindbrain = FindNode(h2vols, "Hindbrain")$volumes,
  #          Cerebellum = FindNode(h2vols, "Cerebellum")$volumes) %>%
  #   pivot_longer(cols=-Strain) %>%
  #   group_by(name) %>%
  #   summarise(h2=h22(value, Strain)) %>%
  #   mutate(type=ifelse(name %in% c("Cerebral cortex", "Hindbrain", "Cerebellum"), "Volume", "Behaviour")) %>% 
  #   mutate(name=fct_reorder(name, h2)) %>%
  #   arrange(h2) %>% 
  #   ggplot() + 
  #   aes(x=name, y=h2, fill=type) + 
  #   geom_col() + 
  #   scale_fill_brewer(palette = "Set1") + 
  #   xlab("") + 
  #   theme_classic() + 
  #   ggtitle("C) comparative heritability of brain and behaviour") + 
  #   theme(axis.text.x = element_text(angle=45, hjust=1))
  
  tmp3 <- ToDataFrameTree(h2vols, "normh2", "name")
  tmp2 <- brainbehaviour %>% 
    #select(Strain, `Accelerated Rotarod:  Motor Performance Improvements`:`Surface righting`) %>% 
    #select(Strain, `Accelerated Rotarod:  Motor Performance Improvements`:`Skilled Reaching: Learning Rate`) %>%
    #select(-Sex.y) %>% 
    dplyr::select(Strain, `Accelerated Rotarod:  Motor Performance Improvements`:Velocity) %>% 
    dplyr::select(-matches("Weight"), -matches("Male/Female")) %>%
    mutate(across(-Strain, as.numeric)) %>% 
    pivot_longer(cols=-Strain) %>% 
    group_by(name) %>%
    summarise(h2=h22(value, Strain))
  tmp1 <- ToDataFrameTree(h2vols, "absh2", "name")
  pH22 <- rbind(tmp1 %>% mutate(h2=absh2, type="Abs. Volume") %>% dplyr::select(h2, name, type), 
                tmp3 %>% mutate(h2=normh2, type="Norm. Volume") %>% dplyr::select(h2, name, type), 
                tmp2 %>% mutate(type="Behaviour")) %>% 
    mutate(type=fct_relevel(type, "Abs. Volume", "Norm. Volume", "Behaviour")) %>% 
    ggplot() + 
    aes(x=h2, fill=type) + 
    geom_density(alpha=0.5) + 
    scale_fill_brewer(palette="Paired") + 
    theme_classic() + 
    ggtitle("C) comparative heritability of brain and behaviour")

  pBV + (pSlices / pH22 + plot_layout(heights=c(7,3))) + plot_layout(widths=c(4,6))
}