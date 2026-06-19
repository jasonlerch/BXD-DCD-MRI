combineMRIbehaviour <- function(newcsv, oldcsv, motorlearning, motorcontrol) {
  gfnew <- read.csv(newcsv) %>% mutate(bn=basename(Filenames),
                                       bn=str_replace(bn, '(.+)_distcorr.+', '\\1'))
  gfold <- read.csv(oldcsv) %>% mutate(bn=basename(Filenames),
                                       bn=str_replace(bn, '(.+)_distcorr.+', '\\1'),
                                       smallID=str_split_i(mouseID, "-", 2),
                                       smallID=as.numeric(smallID))
  mole <- read_excel(motorlearning, skip=4)
  moco <- read_excel(motorcontrol, sheet="Formatted Data Sheet", skip=4)
  
  mri <- gfnew %>%
    left_join(gfold, by="bn", suffix=c(".new", ".old")) %>%
    mutate(motor=ifelse(is.na(motor), "no", motor),
           MRI=TRUE,
           Genotype = str_replace(Genotype, "B6", "C57BL_6J"),
           Genotype = str_replace(Genotype, "DBA", "DBA_2J"),
           Strain = Genotype)

  behav <- full_join(mole, moco, by=c("ID", "Strain")) %>%
    mutate(Strain = str_replace_all(Strain, "0", ""),
           smallID=str_split_i(ID, "_", 1),
           smallID=as.numeric(smallID),
           MRI=FALSE)
  
  full_join(mri, behav, by=c("Strain", "smallID"), na_matches="never") %>%
    mutate(MRI = ifelse(is.na(MRI.x), FALSE, TRUE))
}

combineMRIbehaviour2025 <- function(newcsv, oldcsv, motorheritability, motorcoordination) {
  gfnew <- read.csv(newcsv) %>% mutate(bn=basename(Filenames),
                                       bn=str_replace(bn, '(.+)_distcorr.+', '\\1'))
  gfold <- read.csv(oldcsv) %>% mutate(bn=basename(Filenames),
                                       bn=str_replace(bn, '(.+)_distcorr.+', '\\1'),
                                       smallID=str_split_i(mouseID, "-", 2),
                                       smallID=as.numeric(smallID))
  mri <- gfnew %>%
    left_join(gfold, by="bn", suffix=c(".new", ".old")) %>%
    mutate(motor=ifelse(is.na(motor), "no", motor),
           MRI=TRUE,
           Genotype = str_replace(Genotype, "B6", "C57BL_6J"),
           Genotype = str_replace(Genotype, "DBA", "DBA_2J"),
           Strain = Genotype)
  
  motorherit <- read_excel(motorheritability, range="A5:N293") %>% 
    mutate(smallID= as.numeric(str_split_i(ID, "_", 1)), 
           Strain = str_replace_all(Strain, "BXD0{1,2}(\\d{1,2})", "BXD\\1"))
  
  gait <- read_excel(motorcoordination, sheet="Gait", range="A1:J251")
  rotarod <- read_excel(motorcoordination, sheet="Rotarod", range="A1:E292")
  openfield <- read_excel(motorcoordination, sheet="Openfield", range="A1:I295")
  mcoord <- reduce(list(gait, rotarod, openfield), inner_join, by=c("ID"="Strain"))#inner_join(gait, rotarod, by=c("ID"="Strain"))
  
  motor <- full_join(motorherit, mcoord, by="ID")
  
  mri2 <- left_join(mri, motor, by=c("Strain", "smallID"))
  
  # make numeric columns numeric and identify which group they were in
  mri2 <- mri2 %>% 
    mutate(across(`Accelerated Rotarod:  Motor Performance Improvements`:Velocity & 
                    -matches("Weight") & -matches("Male"), as.numeric)) %>%
    mutate(group = 
             case_when(motor == "no" ~ "no motor training", 
                       motor != "no" & is.na(`Complex Wheel:  Latency to Fall`) & 
                         is.na(`Skilled Reaching: First Time Success`) ~ "non-specific motor training", 
                       !is.na(`Complex Wheel:  Latency to Fall`) ~ "complex wheel", 
                       !is.na(`Skilled Reaching: First Time Success`) ~ "skilled reaching"),
           learners = ifelse(Strain %in% c("BXD27", "BXD28", "BXD75", "BXD15", "BXD86"), "poor", "good"))
  
  return(mri2)
}  