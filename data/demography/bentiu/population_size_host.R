#' Get data from UNWPP
fread("./data/demography/unwpp_population_data.csv") %>%
  melt(id.vars = c("country", "year"), variable.name = "age") %>%
  .[, value := value * 1000] %>%
  .[, age := as.numeric(as.character(age))] %>%
  .[country == "South Sudan" & year == 2015] %>%
  cbind(setAgeBreaks(c(1:100), maxage = set_units(120, "years")), .) %>%
  .[, -c("country", "year", "age")] %>%
  #assume host population data is 70K
  .[, value := value * 7e4/sum(value)] %>%
  .[]
