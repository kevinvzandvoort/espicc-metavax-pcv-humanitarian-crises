#' Get data from UNWPP
fread("./data/demography/unwpp_population_data.csv") %>%
  melt(id.vars = c("country", "year"), variable.name = "age") %>%
  .[, value := value * 1000] %>%
  .[, age := as.numeric(as.character(age))] %>%
  .[country == "Nigeria" & year == 2019] %>%
  cbind(setAgeBreaks(c(1:100), maxage = set_units(120, "years")), .) %>%
  .[, -c("country", "year", "age")] %>%
  #assume IDP population size is 300K
  .[, value := value * 3e5/sum(value)] %>%
  .[]