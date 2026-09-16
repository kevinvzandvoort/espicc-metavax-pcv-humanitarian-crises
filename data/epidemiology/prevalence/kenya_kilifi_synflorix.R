fread("./data/epidemiology/prevalence/carriage_kenya_kilifi.csv") %>%
  cbind(setAgeBreaks(c(0,1,6,15,20,50))) %>%
  setorder(name) %>%
  .[, S := N-VT-NVT] %>%
  .[, -"N"] %>%
  .[, c("name", "from", "to", "S", "VT", "NVT")]
