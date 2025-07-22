fread("./data/epidemiology/prevalence/digaale_pneumosil_alternative_data.csv") %>%
  setNames(c("agegrp", "S", "NVT", "NVT2", "VT", "VT2", "B")) %>%
  .[, agegrp := factor(agegrp, c("<2", "2-5", "6-14", "15-29", "30-49", "50+"))] %>%
  setorder(agegrp) %>%
  .[, c("S", "VT", "NVT", "VT2", "B", "NVT2")] %>%
  cbind(setAgeBreaks(c(0, 2, 6, 15, 30, 50)))
