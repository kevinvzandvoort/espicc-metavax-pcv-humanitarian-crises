pacman::p_load(data.table, socialmixr, janitor)
contact_data = fread("./data/demography/kilifi/contacts.csv")

contact_data[, part_age_group := factor(
  age_class_part,
  c(1:6),
  c("[0y, 1y)", "[1y, 6y)", "[6y, 15y)", "[15y, 20y)", "[20y, 50y)", "[50y, 120y)")
  )]

contact_data[, cont_age_group := factor(
  age_class_cont,
  c("Infant", "Pre-School", "Primary", "Secondary", "Adult", "Older"),
  c("[0y, 1y)", "[1y, 6y)", "[6y, 15y)", "[15y, 20y)", "[20y, 50y)", "[50y, 120y)")
  )]

contact_data_intra = contact_data[share_hh == "Yes", c("csid", "part_age_group", "cont_age_group")]
contact_data_extra = contact_data[share_hh == "No", c("csid", "part_age_group", "cont_age_group")]

contact_data_intra = contact_data_intra %>%
  .[, .(contacts = .N), by = c("part_age_group", "cont_age_group")] %>%
  merge(contact_data_intra[, .(total = length(unique(csid))), by="part_age_group"],
        by = "part_age_group") %>%
  .[, value := contacts/total, by = c("part_age_group", "cont_age_group")]

contact_data_extra = contact_data_extra %>%
  .[, .(contacts = .N), by = c("part_age_group", "cont_age_group")] %>%
  merge(contact_data_extra[, .(total = length(unique(csid))), by="part_age_group"],
        by = "part_age_group") %>%
  .[, value := contacts/total, by = c("part_age_group", "cont_age_group")]

contact_age_groups = setAgeBreaks(c(0, 1, 6, 15, 20, 50))

#' adjust for Kenyan population size in 2011
population_data = fread("./data/demography/kilifi/unpopulation_dataportal_kenya2011.csv") %>%
  janitor::clean_names()
population_data[age_end == "null", age_end := "120"]
population_data[, age_end := as.numeric(age_end)]

for(i in 1:nrow(contact_age_groups)){
  population_data[age_start >= contact_age_groups[i, as.numeric(from)] & age_end <= contact_age_groups[i, as.numeric(to)],
                  age_group := contact_age_groups[i, name]]
}
population_data = population_data[, .(value = sum(value)), by="age_group"]

contact_data_intra = contact_data_intra %>%
  merge(population_data, by.x = "part_age_group", by.y = "age_group") %>%
  .[, total := value.x * value.y] %>%
  .[, c("part_age_group", "cont_age_group", "total")]
contact_data_extra = contact_data_extra %>%
  merge(population_data, by.x = "part_age_group", by.y = "age_group") %>%
  .[, total := value.x * value.y] %>%
  .[, c("part_age_group", "cont_age_group", "total")]

contact_data_intra = contact_data_intra %>%
  merge(contact_data_intra,
        by.x = c("part_age_group", "cont_age_group"),
        by.y = c("cont_age_group", "part_age_group")) %>%
  .[, value := (total.x + total.y)/2] %>%
  merge(population_data, by.x = "part_age_group", by.y = "age_group") %>%
  .[, value := value.x/value.y] %>%
  .[, c("part_age_group", "cont_age_group", "value")]

contact_data_extra = contact_data_extra %>%
  merge(contact_data_extra,
        by.x = c("part_age_group", "cont_age_group"),
        by.y = c("cont_age_group", "part_age_group")) %>%
  .[, value := (total.x + total.y)/2] %>%
  merge(population_data, by.x = "part_age_group", by.y = "age_group") %>%
  .[, value := value.x/value.y] %>%
  .[, c("part_age_group", "cont_age_group", "value")]

list(
  contact_matrix_intra = contact_data_intra %>%
    dcast(part_age_group~cont_age_group) %>%
    .[, -"part_age_group"] %>%
    as.matrix(),
  contact_matrix_extra = contact_data_extra %>%
  dcast(part_age_group~cont_age_group) %>%
  .[, -"part_age_group"] %>%
  as.matrix(),
  population = population_data %>%
    merge(contact_age_groups, by.x="age_group", by.y = "name") %>%
    .[, c("from", "value")] %>%
    setNames(c("lower.lage.limit", "population")))
