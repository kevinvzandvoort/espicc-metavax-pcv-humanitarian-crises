#' Age groups that will be used in the model
age_groups_model = c(set_units(c(0:12), "month"),
                     set_units(c(1:5, 6, 10, 15, 20, 30, 40, 50), "years")) %>%
  setAgeBreaks() %>% .[, age := .I]

modelSetup = function(contact_data = list(data$demography$digaale$displaced$contact_data$data,
                                              data$demography$digaale$host$contact_data$data),
                      contacts_extra_host = data$demography$digaale$displaced$contacts_extra_host$data$value,
                      contacts_prop_extra = data$demography$digaale$displaced$contacts_prop_extra$data$value,
                      population_data = list(data$demography$digaale$displaced$population_data$data,
                                                       data$demography$digaale$host$population_data$data),
                      prop_malnourished = c(data$demography$digaale$displaced$malnourished$data$value,
                                            data$demography$digaale$host$malnourished$data$value),
                      migration_rate = data$demography$digaale$displaced$migration_rate$data$value,
                      malnourished_rr_transmission = data$epidemiology$malnourished_RR_transmission$data$value,
                      malnourished_rr_disease = data$epidemiology$malnourished_RR_disease$data$value,
                      clearance_rate = data$epidemiology$clearance_rate$data,
                      vaccination_groups = c("unvaccinated", "pcv1", "pcv2")){
  
  #' increased acquisition rate in those malnourished
  maln_rr_transmission = setAgeBreaks(0) %>%
    .[, value := data$epidemiology$malnourished_RR_transmission$data$value] %>%
    combineAgeBreaks(x = age_groups_model, y = .) %>% .[, value]
  
  #' risk ratio of 1 in all age groups
  null_rr = setAgeBreaks(0) %>%
    .[, value := 1] %>%
    combineAgeBreaks(x = age_groups_model, y = .) %>% .[, value]
  
  #' clearance rates per day
  clearance_rates = clearance_rate
  
  #' The competition parameter reflects the extent to which carrying one or more pneumococcal serotypes of one group (carrying
  #' VT or NVT) protects against acquisition of serotypes in the other group (NVT or VT). I.e. if competition is 20%, those
  #' who carry VT serotypes experience 0.8 times the rate of infection with NVTs compared to someone who is susceptible
  #' NB - is serotype specific in the real world
  pneumo_competition = 0.90 #fit this
  
  #' D - Population specific parameters
  #displaced_population_size_cm = digaale_contact_matrix_agegroups %>% combineAgeBreaks(population_data[[1]], method = "sum")
  displaced_population_size_model = age_groups_model %>% combineAgeBreaks(population_data[[1]], method = "sum")
  displaced_contact_matrix_model = adjustContactMatrixAgeGroups(age_groups_model, contact_data[[1]],
                                                                {if(length(population_data) == 2) setAgeBreaks(c(0, 10, 20, 30, 40, 50)) else setAgeBreaks(seq(0, 75, 5))},
                                                                displaced_population_size_model)
  
  #proportion of extra-household contacts made with host population
  if(length(population_data) == 2){
    populations = c("idp_malnourished", "idp_non_malnourished", "host_malnourished", "host_non_malnourished")
    
    host_population_size_model = age_groups_model %>% combineAgeBreaks(population_data[[2]], method = "sum")
    host_contact_matrix_model = adjustContactMatrixAgeGroups(age_groups_model, contact_data[[2]], setAgeBreaks(seq(0, 75, 5)), host_population_size_model)
    #' ensure the matrices have the same dominant eigenvalue
    host_contact_matrix_model = host_contact_matrix_model*(dominantEigenValue(displaced_contact_matrix_model)/dominantEigenValue(host_contact_matrix_model))
    
    displaced_prop_contacts_host = contacts_extra_host * contacts_prop_extra
    travel_matrix = calculateTravelMatrix(displaced_prop_contacts_host,
                                          displaced_contact_matrix_model, displaced_population_size_model, prop_malnourished[1],
                                          host_contact_matrix_model, host_population_size_model, prop_malnourished[2])
    
    migration_matrix = calculateMigrationMatrix(migration_rate)
  } else {
    populations = c("idp_malnourished", "idp_non_malnourished")
    
    travel_matrix = matrix(c(prop_malnourished, 1-prop_malnourished, prop_malnourished, 1-prop_malnourished),
                           nrow = 2, ncol = 2, dimnames = list("contactee" = c("malnourished", "nourished"),
                                                               "contactor" = c("malnourished", "nourished")))
    migration_matrix = matrix(c(0, 0, 0, 0), nrow = 2, ncol = 2,
                              dimnames = list("host" = c("malnourished", "nourished"),
                                              "migrant" = c("malnourished", "nourished")))
  }
  
  #' These will be used to replace the arms lists in the model_populations
  #' Values repeated in vaccination_strategies
  no_coverage = list(
    list(value = getVaccineCoverage(age_groups_model, set_units(0, "months"), 0),
         time = 0, coverage_to = rep(1, age_groups_model[, .N])))
  ve_carriage = setAgeBreaks(0) %>% .[, value := 0] %>%
    combineAgeBreaks(x = age_groups_model, y = .) %>% .[, value]
  ve_disease = setAgeBreaks(0) %>% .[, value := 0] %>%
    combineAgeBreaks(x = age_groups_model, y = .) %>% .[, value]
  vwaning = setAgeBreaks(0) %>% .[, value := 0] %>%
    combineAgeBreaks(x = age_groups_model, y = .) %>% .[, value]
  
  #' Create model_params object and specify populations to model
  #' - Make sure to transpose the contact matrices so that columns denote contactees and rows contactors
  #'   as this is required for the matrix multiplication when calculating the FOIs in the model
  #'   i.e row i will show average number of contacts made by those of age i in all age groups, which will
  #'   be multiplied with a column vector with prevalence in each contacted age groups and summed to give
  #'   the average number of effective contacts with infectious individuals made by someone aged i (or j in
  #'   contact matrix before transposing)
  #' - We adjust the probability of acquiring pneumococci for those who are malnourished, but not for those
  #'   not malnourished (NULL_RR vector has 1 in all age groups)
  #' - Nb need to refactor start and stop of acq adjustments. Set for negative and very large timesteps now
  #'   to make sure they are always active
  model_params = model_params %>%
    createParameters(age_groups_model,
                     populations = populations,
                     vaccination_groups = vaccination_groups) %>%
    setParameter(level = "global", key = "initial_states",
                 value = c(0.8, 0.1, 0.1, 0) %>%
                   rep(each = age_groups_model[, .N]) %>%
                   c(c(0, 0, 0, 0) %>% rep(each = age_groups_model[, .N]) %>%
                       rep(length(vaccination_groups) - 1)) %>%
                   rep(length(populations))) %>%
    setParameter(level = "global", key = "comp", value = pneumo_competition) %>%
    setParameter(level = "global", key = "clearVT",
                 value = age_groups_model %>% combineAgeBreaks(clearance_rates[st == "VT"]) %>% .[, value]) %>%
    setParameter(level = "global", key = "clearNVT",
                 value = age_groups_model %>% combineAgeBreaks(clearance_rates[st == "NVT"]) %>% .[, value]) %>%
    setParameter(level = "global", key = "travel", value = travel_matrix) %>%
    setParameter(level = "global", key = "migration", value = migration_matrix) %>%
    setParameter(level = "population", population = "idp_malnourished", key = "N",
                 value = displaced_population_size_model[, value] * prop_malnourished[1]) %>%
    setParameter(level = "population", population = "idp_non_malnourished", key = "N",
                 value = displaced_population_size_model[, value] * (1 - prop_malnourished[1])) %>%
    setParameter(level = "population", population = "all", key = "adjust_acq_start", value = -1) %>%
    setParameter(level = "population", population = "all", key = "adjust_acq_stop", value = 1e6) %>%
    setParameter(level = "population", population = c("idp_malnourished"),
                 key = "adjust_acq", value = maln_rr_transmission) %>%
    setParameter(level = "population", population = c("idp_non_malnourished"),
                 key = "adjust_acq", value = null_rr) %>%
    setParameter(level = "population", population = c("idp_malnourished", "idp_non_malnourished"),
                 key = "betaVT", value = displaced_contact_matrix_model %>% t) %>%
    setParameter(level = "population", population = c("idp_malnourished", "idp_non_malnourished"),
                 key = "betaNVT", value = displaced_contact_matrix_model %>% t) %>%
    setParameter(level = "vaccination_group", population = "all", vaccination_group = "all",
                 key = "coverage_r", value = no_coverage) %>%
    setParameter(level = "vaccination_group", population = "all", vaccination_group = "all",
                 key = "coverage_c", value = no_coverage) %>%
    setParameter(level = "vaccination_group", population = "all", vaccination_group = "all",
                 key = "efficacy_transmission", value = rep(0, age_groups_model[, .N])) %>%
    setParameter(level = "vaccination_group", population = "all", vaccination_group = "all",
                 key = "efficacy_disease", value = rep(0, age_groups_model[, .N])) %>%
    setParameter(level = "vaccination_group", population = "all", vaccination_group = c("pcv1", "pcv2"),
                 key = "waning", value = vwaning, strict = FALSE)
  
  #' the host population is excluded in the Bambari scenario
  if(length(population_data) == 2){
    model_params = model_params %>%
      setParameter(level = "population", population = "host_malnourished", key = "N",
                   value = host_population_size_model[, value] * prop_malnourished[2]) %>%
      setParameter(level = "population", population = "host_non_malnourished", key = "N",
                   value = host_population_size_model[, value] * (1 - prop_malnourished[2])) %>%
      setParameter(level = "population", population = c("host_malnourished"),
                   key = "adjust_acq", value = maln_rr_transmission) %>%
      setParameter(level = "population", population = c("host_non_malnourished"),
                   key = "adjust_acq", value = null_rr) %>%
      setParameter(level = "population", population = c("host_malnourished", "host_non_malnourished"),
                   key = "betaVT", value = displaced_contact_matrix_model %>% t) %>%
      setParameter(level = "population", population = c("host_malnourished", "host_non_malnourished"),
                   key = "betaNVT", value = displaced_contact_matrix_model %>% t)
  }
  
  #' use difference equations
  model_params = model_params %>%
    setParameter(key = "model_solver_type", value = "DIFF") %>%
    setParameter(key = "model_solver_difference", value = TRUE) %>%
    setParameter(key = "solver_difference_delta_t", value = 1)
  
  #' add any additional values needed in updateParameters
  model_params$malnourished_RR_disease = data$epidemiology$malnourished_RR_disease$data$value
  
  if(length(populations) == 4){
    host_contact_matrix_model = host_contact_matrix_model
    model_params$cm_unadjusted = list(displaced = displaced_contact_matrix_model,
                                      host = host_contact_matrix_model)
  } else {
    model_params$cm_unadjusted = list(displaced = displaced_contact_matrix_model)  
  }
  
  return(model_params)
}

#' Parameters to fit
priors = rbindlist(list(
  data.table(variable = "beta_1",
             min = 1e-20, max = 1,
             plotgroup = "beta_values",
             density = function(x, uselog=TRUE) dbeta(x, shape1 = 0.1, shape2 = 10, log=uselog),
             sampler = function(n) rbeta(n, shape1 = 0.1, shape2 = 10)),
  data.table(variable = "beta_2", 
             min = 1e-20, max = 1,
             plotgroup = "beta_values",
             density = function(x, uselog=TRUE) dbeta(x, shape1 = 0.1, shape2 = 10, log=uselog),
             sampler = function(n) rbeta(n, shape1 = 0.1, shape2 = 10)),
  data.table(variable = "beta_3",
             min = 1e-20, max = 1,
             plotgroup = "beta_values",
             density = function(x, uselog=TRUE) dbeta(x, shape1 = 0.1, shape2 = 10, log=uselog),
             sampler = function(n) rbeta(n, shape1 = 0.1, shape2 = 10)),
  data.table(variable = "beta_VT_NVT_u5",
             min = 0.2, max = 5,
             plotgroup = "beta_ratio_values",
             density = function(x, uselog=TRUE) dlnorm(x, meanlog = log(1^2 / sqrt(0.2^2 + 1^2)), sdlog = sqrt(log(1 + 0.2^2 / 1^2)), log=uselog),
             sampler = function(n) rlnorm(n, meanlog = log(1^2 / sqrt(0.2^2 + 1^2)), sdlog = sqrt(log(1 + 0.2^2 / 1^2)))),
  data.table(variable = "beta_VT_NVT_o5",
             min = 0.2, max = 5,
             plotgroup = "beta_ratio_values",
             density = function(x, uselog=TRUE) dlnorm(x, meanlog = log(1^2 / sqrt(0.2^2 + 1^2)), sdlog = sqrt(log(1 + 0.2^2 / 1^2)), log=uselog),
             sampler = function(n) rlnorm(n, meanlog = log(1^2 / sqrt(0.2^2 + 1^2)), sdlog = sqrt(log(1 + 0.2^2 / 1^2)))),
  data.table(variable = "competition",
             min = 1e-20, max = 1,
             plotgroup = "competition",
             density = function(x, uselog=TRUE) dbeta(x, shape1 = 4.53, shape2 = 17.14, log=uselog),
             sampler = function(n) rbeta(n, shape1 = 4.53, shape2 = 17.14))))

priors[, variable := factor(variable, variable)]

#' Create prior for use with BayesianTools
prior = createBTPrior(priors)

additional_posteriors = rbindlist(list(
  data.table(variable = "efficacy_carriage", sampler = function(n) rbeta(n, shape1 = 47.29984, shape2 = 47.29984)),
  data.table(variable = "efficacy_total", sampler = function(n) rbeta(n, shape1 = 12, shape2 = 3)),
  data.table(variable = "efficacy_duration", sampler = function(n) rlnorm(n, meanlog = log(6), sdlog = 0.25)),
  data.table(variable = "ccr_VT_u1", sampler = function(n) rlnorm(n, meanlog = -7.287789, sdlog = 0.2440571)),
  data.table(variable = "ccr_VT_1_5", sampler = function(n) rlnorm(n, meanlog = -8.342471, sdlog = 0.1491026)),
  data.table(variable = "ccr_VT_o5", sampler = function(n) rlnorm(n, meanlog = -10.827291, sdlog = 0.4109898)),
  data.table(variable = "ccr_NVT_u1", sampler = function(n) rlnorm(n, meanlog = -8.811200, sdlog = 0.2083088)),
  data.table(variable = "ccr_NVT_1_5", sampler = function(n) rlnorm(n, meanlog = -10.718420, sdlog = 0.1937705)),
  data.table(variable = "ccr_NVT_o5", sampler = function(n) rlnorm(n, meanlog = -12.618624, sdlog = 0.2478390))))

#' Process additional posteriors to combine with BayesianTools prior
additional_posteriors = createBTAdditionalPosterior(additional_posteriors)

#showPriorsTable(priors)[, type := "fitted"] %>%
#  rbind(showPriorsTable(additional_posteriors$setup)[, type := "fixed"]) %>%
#  fwrite("prior_values_table.csv")

updateDepParameters = function(params){
  #' make sure that efficacy_total can't be lower than efficacy_carriage
  if(all(c("efficacy_total", "efficacy_carriage") %in% names(params))){
    if(params["efficacy_total"] < params["efficacy_carriage"]){
      params["efficacy_carriage"] = params["efficacy_total"]
    } 
  }
  
  return(params)
}

#' Function that will be used to replace the model parameters with new values
updateParameters = function(params, model_params){
  params = updateDepParameters(params)
  
  age_groups_model = model_params$global_settings$age_groups_model
  
  #' Calculate new contact matrices
  betaVT = set_units(c(0, 5, 15), "year") %>% setAgeBreaks() %>%
    .[, value := c(params["beta_1"]*(params["beta_VT_NVT_u5"]),
                   params["beta_2"]*(params["beta_VT_NVT_o5"]),
                   params["beta_3"]*(params["beta_VT_NVT_o5"]))] %>%
    combineAgeBreaks(x = age_groups_model, y = .) %>% .[, value]
  
  betaNVT = set_units(c(0, 5, 15), "year") %>% setAgeBreaks() %>%
    .[, value := c(params["beta_1"],
                   params["beta_2"],
                   params["beta_3"])] %>%
    combineAgeBreaks(x = age_groups_model, y = .) %>% .[, value]
  
  new_matrix_displaced_VT = sweep(model_params$cm_unadjusted$displaced, 2, betaVT, "*") %>% t
  new_matrix_displaced_NVT = sweep(model_params$cm_unadjusted$displaced, 2, betaNVT, "*") %>% t
  
  model_params = model_params %>%
    setParameter(level = "population", population = c("idp_malnourished", "idp_non_malnourished"),
                 key = "betaVT", value = new_matrix_displaced_VT) %>%
    setParameter(level = "population", population = c("idp_malnourished", "idp_non_malnourished"),
                 key = "betaNVT", value = new_matrix_displaced_NVT)
  
  #' add host-population contact matrices, if used
  if(length(model_params$populations) == 4){
    new_matrix_host_VT = sweep(model_params$cm_unadjusted$host, 2, betaVT, "*") %>% t
    new_matrix_host_NVT = sweep(model_params$cm_unadjusted$host, 2, betaNVT, "*") %>% t
    
    model_params = model_params %>%
      setParameter(level = "population", population = c("host_malnourished", "host_non_malnourished"),
                   key = "betaVT", value = new_matrix_host_VT) %>%
      setParameter(level = "population", population = c("host_malnourished", "host_non_malnourished"),
                   key = "betaNVT", value = new_matrix_host_NVT)
  }
  
  #' add new competition value
  model_params = model_params %>%
    setParameter(level = "global", key = "comp", value = params["competition"])
  
  #' add new vaccine values
  #' - assume VE is 50% lower in infants only receiving a single dose of PCV
  ve_carriage = set_units(c(0, 1), "year") %>% setAgeBreaks() %>%
    .[, value := c(params["efficacy_carriage"], params["efficacy_carriage"])] %>%
    combineAgeBreaks(x = age_groups_model, y = .) %>% .[, value]
  ve_disease = set_units(c(0, 1), "year") %>% setAgeBreaks() %>%
    .[, value := c(params["efficacy_total"], params["efficacy_total"])] %>%
    combineAgeBreaks(x = age_groups_model, y = .) %>% .[, value]
  ve_disease = 1 - (1 - ve_disease)/(1 - ve_carriage)
  ve_waning = set_units(c(0), "year") %>% setAgeBreaks() %>%
    .[, value := 1/(params["efficacy_duration"] * 365)] %>%
    combineAgeBreaks(x = age_groups_model, y = .) %>% .[, value]
  
  model_params = model_params %>%
    setParameter(level = "vaccination_group", population = "all", vaccination_group = "pcv1",
                 key = "efficacy_transmission", value = ve_carriage * 0.5, strict = FALSE) %>%
    setParameter(level = "vaccination_group", population = "all", vaccination_group = "pcv2",
                 key = "efficacy_transmission", value = ve_carriage, strict = FALSE) %>%
    setParameter(level = "vaccination_group", population = "all", vaccination_group = "pcv1",
                 key = "efficacy_disease", value = ve_disease * 0.5, strict = FALSE) %>%
    setParameter(level = "vaccination_group", population = "all", vaccination_group = "pcv2",
                 key = "efficacy_disease", value = ve_disease, strict = FALSE) %>%
    setParameter(level = "vaccination_group", population = "all", vaccination_group = c("pcv1", "pcv2"),
                 key = "waning", value = ve_waning, strict = FALSE)
  
  #' add new case-carrier ratios
  ccrNVT = set_units(c(0, 1, 5), "year") %>% setAgeBreaks() %>%
    .[, value := c(params["ccr_NVT_u1"],
                   params["ccr_NVT_1_5"],
                   params["ccr_NVT_o5"])] %>%
    combineAgeBreaks(x = age_groups_model, y = .) %>% .[, value]
  ccrVT = set_units(c(0, 1, 5), "year") %>% setAgeBreaks() %>%
    .[, value := c(params["ccr_VT_u1"],
                   params["ccr_VT_1_5"],
                   params["ccr_VT_o5"])] %>%
    combineAgeBreaks(x = age_groups_model, y = .) %>% .[, value]
  model_params = model_params %>%
    setParameter(level = "population", population = "idp_non_malnourished",
                 key = "ccrNVT", value = ccrNVT) %>%
    setParameter(level = "population", population = "idp_non_malnourished",
                 key = "ccrVT", value = ccrVT) %>%
    setParameter(level = "population", population = "idp_malnourished",
                 key = "ccrNVT", value = ccrNVT * model_params$malnourished_RR_disease) %>%
    setParameter(level = "population", population = "idp_malnourished",
                 key = "ccrVT", value = ccrVT * model_params$malnourished_RR_disease)
  
  #' add host-population CCRs, if used
  if(length(model_params$populations) == 4){
    model_params = model_params %>%
      setParameter(level = "population", population = "host_non_malnourished",
                   key = "ccrNVT", value = ccrNVT) %>%
      setParameter(level = "population", population = "host_non_malnourished",
                   key = "ccrVT", value = ccrVT) %>%
      setParameter(level = "population", population = "host_malnourished",
                   key = "ccrNVT", value = ccrNVT * model_params$malnourished_RR_disease) %>%
      setParameter(level = "population", population = "host_malnourished",
                   key = "ccrVT", value = ccrVT * model_params$malnourished_RR_disease)
  }
  
  return(model_params)
}
