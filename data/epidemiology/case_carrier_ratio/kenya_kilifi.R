fread("./data/epidemiology/case_carrier_ratio/flasche_kenya_ccr.csv") %>%
  merge(setAgeBreaks(c(0, 1, 6, 15, 21, 50)) %>%
        .[, age_group := c("<1y", "1-5y", "6-14y", "15-20y", "21-49y", "50+")],
      by="age_group") %>%
  .[, c("iter", "st", "from", "to", "name", "value")] %>%
  setorder(iter, st, name) %>%
  merge(setAgeBreaks(c(0, 1, 6, 15, 21, 50)) %>%
          .[, age_group := c("<1y", "1-5y", "6-14y", "15-20y", "21-49y", "50+")],
        by="age_group") %>%
  .[, c("iter", "st", "from", "to", "name", "value")] %>%
  setorder(iter, st, name)