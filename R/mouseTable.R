makeMouseTable <- function(brainbehaviourNoOut) {
  brainbehaviourNoOut %>% 
    mutate(Genotype = str_replace_all(Genotype, "_", "/")) %>% 
    mutate(group=fct_relevel(group, "no motor training", "non-specific motor training")) %>% 
    group_by(Genotype, group) %>% 
    summarise(n=n()) %>% 
    pivot_wider(values_from = n, names_from = group) %>% 
    ungroup() %>% gt(rowname_col = "Genotype") %>% 
    sub_missing() %>% 
    tab_header("Mouse numbers") %>% 
    summary_columns(fns=list("sum"), new_col_names = "Total")  %>% 
    grand_summary_rows(fns=list(list(id="sum", label="Total") ~ sum(., na.rm = T)))
}

mouseTableCategory <- function(brainbehaviourNoOut) {
  brainbehaviourNoOut %>% 
    select(mouseID, 
           Genotype, 
           `Accelerated Rotarod:  Motor Performance Improvements`:Velocity, 
           -`Male/Female.x`, 
           -`Male/Female.y`, 
           -`Male/Female`) %>% 
    pivot_longer(`Accelerated Rotarod:  Motor Performance Improvements`:Velocity) %>% 
    mutate(bt=assignBehaviourTest(name)) %>% 
    group_by(mouseID, bt) %>% 
    summarise(n=sum(any(!is.na(value))), 
              Genotype=Genotype[1]) %>% 
    group_by(Genotype, bt) %>% summarise(n=sum(n)) %>% 
    filter(!is.na(bt)) %>% 
    pivot_wider(names_from = bt, values_from = n) %>%
    ungroup() %>%
    mutate(Genotype = str_replace_all(Genotype, "_", "/")) %>% 
    left_join(mouseTableBehaviour(brainbehaviourNoOut)) %>%
    gt(rowname_col = "Genotype") %>% 
    sub_missing() %>% 
    cols_move_to_start(c("Total", "control", "motor testing")) %>%
    tab_spanner("Overall numbers", level=2, c("Total", "control", "motor testing")) %>%
    tab_spanner("Motor Learning", c("Accelerated Rotarod", "Complex Wheel", "Horizontal Rung", "Skilled Reaching")) %>%
    tab_spanner("Generalized test of motor function", c("Gait", "Open Field", "Rotarod")) %>%
    tab_spanner("Motor test numbers", c("Accelerated Rotarod", "Complex Wheel", "Gait", "Horizontal Rung", "Open Field", "Rotarod", "Skilled Reaching")) %>%
    tab_header("Mouse numbers") %>% 
    #summary_columns(fns=list("sum"), new_col_names = "Total")  %>% 
    grand_summary_rows(fns=list(list(id="sum", label="Total") ~ sum(., na.rm = T)))
}

mouseTableBehaviour <- function(brainbehaviourNoOut) {
  brainbehaviourNoOut %>% 
    mutate(Genotype = str_replace_all(Genotype, "_", "/")) %>% 
    mutate(group=ifelse(group == "no motor training", "control", "motor testing")) %>% 
    group_by(Genotype, group) %>% 
    summarise(n=n()) %>%
    pivot_wider(values_from = n, names_from = group) %>% 
    mutate(Total=sum(c(control, `motor testing`), na.rm=T)) %>%
    ungroup()
}