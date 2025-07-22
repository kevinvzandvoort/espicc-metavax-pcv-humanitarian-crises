pacman::p_load(magrittr, socialmixr, data.table)

if(!file.exists("./data/demography/digaale/digaale_contact_data.RDS")){
  #' Get data from Zenodo
  socialmixr::get_survey("https://zenodo.org/record/5226281#.YR-TzlvTVH6") %>%
    saveRDS("./data/demography/digaale/digaale_contact_data.RDS") 
}

#' The estimated population size in Digaale (for provided age groups)
#'  can manually be downloaded
if(!file.exists("./data/demography/digaale/digaale_survey_population.RDS")){
  digaale_survey_population =
    data.table::fread("https://zenodo.org/record/5226281/files/espicc_somaliland_digaale_survey_population.csv")
  saveRDS(digaale_survey_population, "./data/demography/digaale/digaale_survey_population.RDS")
}

#' Note that weekends fall on Fridays and Saturdays in Somaliland.
#' - The dayofweek variable provided in the dataset has been kept
#'   consistent with R defaults (0: Sunday to 6: Saturday)
readRDS("./data/demography/digaale/digaale_contact_data.RDS") %>%
  (function(x){
    x$participants[, c("dayofweek", "dayofweek_name", "weekend")] %>%
      unique %>%
      setorder(dayofweek) %>%
      #' socialmixr currently assumes the weekend to fall on dayofweek
      #'  6 (Saturday) and 0 (Sunday)
      #' - dayofweek can be manually edited so that Fridays and Saturdays
      #'   are taken as the weekend, if you wish to weight contacts by
      #'   weekday
      .[, dayofweek := ifelse(dayofweek == 6, 0, dayofweek + 1)]
    
    #' separate contacts in intra- and extra-household contacts
    x$contacts %<>%
      .[, household := ifelse(contact_relationship == "Household member", "intra", "extra")]    
    
    digaale_survey_population = readRDS("./data/demography/digaale/digaale_survey_population.RDS") %>%
      .[lower.age.limit >= 50, lower.age.limit := 50] %>%
      .[, .(population=sum(population)), by="lower.age.limit"]
    
    #' The contact matrix can then be constructed as follows
    #' - The provided survey_population can be used to construct a
    #'   population representative matrix for Digaale IDP camp
    #' - As the sample is not self-weighing (oversampling of young
    #'   age groups), it is recommended to apply the survey_weight
    #'   as weights
    #' Note socialmixr's contact matrices show contactors in rows
    #'  and contactees in columns
    digaale_contact_matrix_intra = x %>%
      socialmixr::contact_matrix(survey.pop = digaale_survey_population,
                                 age.limits = digaale_survey_population$lower.age.limit,
                                 symmetric = TRUE, weights = "survey_weight", weigh.dayofweek = TRUE,
                                 filter = list("household" = "intra"))
    
    digaale_contact_matrix_extra = x %>%
      socialmixr::contact_matrix(survey.pop = digaale_survey_population,
                                 age.limits = digaale_survey_population$lower.age.limit,
                                 symmetric = TRUE, weights = "survey_weight", weigh.dayofweek = TRUE,
                                 filter = list("household" = "extra"))
    
    list(contact_matrix_intra = digaale_contact_matrix_intra$matrix,
         contact_matrix_extra = digaale_contact_matrix_extra$matrix,
         population = digaale_survey_population)
    
  })
