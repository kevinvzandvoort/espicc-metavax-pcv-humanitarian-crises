pacman::p_load(magrittr, data.table)

if(!file.exists("./data/digaale_population_size.RDS")){
  x = tempfile(fileext = ".RDS")
  download.file("https://github.com/kevinvzandvoort/espicc-somaliland-digaale-survey-2019/blob/main/data/household_data_members.RDS?raw=true", x, mode=ifelse(Sys.info()["sysname"] == "Windows", "wb", "w"))
  saveRDS(readRDS(x), "./data/digaale_population_size.RDS")
  unlink(x)
}

#' Re-allocate some age groups
readRDS("./data/digaale_population_size.RDS") %>%
  .[, age_group2 := age_group_sample] %>%
  .[infant == TRUE, age_group2 := "<1"] %>%
  .[infant == FALSE & age_group_sample == "<2", age_group2 := "1"] %>%
  .[age_group_60 == "0-9" & age_group_sample == "6-14", age_group2 := "6-9"] %>%
  .[age_group_60 == "10-19" & age_group_sample == "6-14", age_group2 := "10-14"] %>%
  .[age_group_60 == "10-19" & age_group_sample == "15-29", age_group2 := "15-19"] %>%
  .[age_group_60 == "20-29" & age_group_sample == "15-29", age_group2 := "20-29"] %>%
  .[age_group_60 == "30-39", age_group2 := "30-39"] %>%
  .[age_group_60 == "40-49", age_group2 := "40-49"] %>%
  .[age_group_80 == "50-59", age_group2 := "50-59"] %>%
  .[age_group_80 == "60-69", age_group2 := "60-69"] %>%
  .[age_group_80 == "70-79", age_group2 := "70-79"] %>%
  .[age_group_80 == "80+", age_group2 := "80+"] %>%
  #Assume annual population <2 is the same as annual 2-5
  .[age_group2 %in% c("<1", "1", "2-5"), age_group2 := "<5"] %>%
  .[, age_group2 := factor(age_group2, c("<5", "6-9", "10-14", "15-19", "20-29", "30-39",
                                                     "40-49", "50-59", "60-69", "70-79", "80+"),
                           setAgeBreaks(c(0, 6, 10, 15, 20, 30, 40, 50, 60, 70, 80))[, name])] %>%
  .[, .(N = .N), by="age_group2"] %>%
  merge(setAgeBreaks(c(0, 6, 10, 15, 20, 30, 40, 50, 60, 70, 80)), ., by.x="name", by.y="age_group2") %>%
  #adjust observed estimates for estimated ratio inhabited/participated shelters
  .[, value := N * 715/489] %>% .[, -"N"]
