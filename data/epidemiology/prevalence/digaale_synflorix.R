fread("./data/epidemiology/prevalence/digaale_data_for_model.csv") %>%
  .[, .N, by=c("pcv10_type", "agegrp")] %>%
  dcast(agegrp ~ pcv10_type, value.var="N") %>%
  .[, agegrp := factor(agegrp, c("<2", "2-5", "6-14", "15-29", "30-49", "50+"))] %>%
  setorder(agegrp) %>%
  setNames(c("agegrp", "NVT", "S", "VT", "B")) %>%
  .[, c("S", "VT", "NVT", "B")] %>%
  cbind(setAgeBreaks(c(0, 2, 6, 15, 30, 50)))