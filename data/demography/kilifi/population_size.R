#' Get data from UNWPP
fread("./data/demography/kilifi/unpopulation_dataportal_kenya2011.csv") %>%
  .[, age := as.numeric(as.character(AgeStart))] %>%
  cbind(setAgeBreaks(c(1:100), maxage = set_units(120, "years")), .) %>%
  .[, value := Value] %>%
  .[, c("name", "from", "to", "value")]
