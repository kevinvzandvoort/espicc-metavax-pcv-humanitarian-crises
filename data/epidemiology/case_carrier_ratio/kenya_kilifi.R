fread("./data/epidemiology/case_carrier_ratio/ccr_kenya_param_posteriors.txt") %>% .[-1] %>%
  .[, c("CCR VT <1y", "CCR VT 1-5y", "CCR VT 6-14y", "CCR VT 15-20y", "CCR VT 21-49y", "CCR VT 50+",
        "CCR NVT <1y", "CCR NVT 1-5y", "CCR NVT 6-14y", "CCR NVT 15-20y", "CCR NVT 21-49y", "CCR NVT 50+")] %>%
  .[, iter := 1:.N] %>%
  melt(id.vars = "iter") %>%
  .[, st := ifelse(grepl("NVT", variable), "NVT", "VT")] %>%
  .[, age_group := strsplit(as.character(variable), split = " ") %>% sapply("[[", 3)] %>%
  merge(setAgeBreaks(c(0, 1, 6, 15, 21, 50)) %>%
          .[, age_group := c("<1y", "1-5y", "6-14y", "15-20y", "21-49y", "50+")],
        by="age_group") %>%
  .[, c("iter", "st", "from", "to", "name", "value")] %>%
  setorder(iter, st, name)