#' Calculate the travel matrix between the populations
calculateTravelMatrix = function(contacts_prop_idp_host, digaale_contact_matrix_model, digaale_population_size_model, digaale_malnourished,
                                 hargeisa_contact_matrix_model, hargeisa_population_data_model, hargeisa_malnourished){
  #' proportion of contacts in host population that are made with IDP population
  contacts_prop_host_idp =
    sum(digaale_contact_matrix_model %*% diag(digaale_population_size_model[, value]) * contacts_prop_idp_host)/
    sum(hargeisa_contact_matrix_model %*% diag(hargeisa_population_data_model[, value]))
  
  travel = travelm = matrix(0, ncol=4, nrow=4)
  travel[, 1] = c((1-contacts_prop_idp_host), (1-contacts_prop_idp_host), contacts_prop_idp_host, contacts_prop_idp_host)
  travel[, 2] = c((1-contacts_prop_idp_host), (1-contacts_prop_idp_host), contacts_prop_idp_host, contacts_prop_idp_host)
  travel[, 3] = c(contacts_prop_host_idp, contacts_prop_host_idp, (1-contacts_prop_host_idp), (1-contacts_prop_host_idp))
  travel[, 4] = c(contacts_prop_host_idp, contacts_prop_host_idp, (1-contacts_prop_host_idp), (1-contacts_prop_host_idp))
  travelm[1,] = c(digaale_malnourished, digaale_malnourished, digaale_malnourished, digaale_malnourished)
  travelm[2,] = c(1-digaale_malnourished, 1-digaale_malnourished, 1-digaale_malnourished, 1-digaale_malnourished)
  travelm[3,] = c(hargeisa_malnourished, hargeisa_malnourished, hargeisa_malnourished, hargeisa_malnourished)
  travelm[4,] = c(1-hargeisa_malnourished, 1-hargeisa_malnourished, 1-hargeisa_malnourished, 1-hargeisa_malnourished)
  
  travel = travel * travelm
  
  dimnames(travel) = list(contactee=c("Digaale/malnourished", "Digaale/nourished", "Hargeisa/malnourished", "Hargeisa/nourished"),
                          contactor=c("Digaale/malnourished", "Digaale/nourished", "Hargeisa/malnourished", "Hargeisa/nourished"))
  
  return(travel)
}

calculateMigrationMatrix = function(migration_rate){
  #' Where a z is provided (-1), the model will automatically calculate the associated migration rate in order to
  #' keep population sizes in Digaale constant
  m = migration_rate
  z = -1
  
  migration_matrix =
    matrix(c(0, 0, z, 0,
             0, 0, 0, z,
             m, 0, 0, 0,
             0, m, 0, 0),
           byrow = TRUE, nrow = 4,
           dimnames = list(host=c("Digaale/malnourished", "Digaale/nourished", "Hargeisa/malnourished", "Hargeisa/nourished"),
                           migrant=c("Digaale/malnourished", "Digaale/nourished", "Hargeisa/malnourished", "Hargeisa/nourished")))
  
  return(migration_matrix)
}

prepareDataForPosterior = function(res_prevacc, data, stratify_cluster = FALSE){
  bycols = c("run", "compartment", "name.y")
  if(stratify_cluster) bycols = c(bycols, "population")
  
  age_groups_model %>% matchingAgeBreaks(data) %>% .[, age := 1:.N] %>% merge(res_prevacc, by="age") %>%
    .[, .(value = sum(value), N = sum(N), prev = sum(value)/sum(N)), by=bycols] %>%
    (function(x){
      if(stratify_cluster) x %>% dcast(run + cluster + name.y ~ compartment, value.var = "prev")
      else x %>% dcast(run + name.y ~ compartment, value.var = "prev")
    }) %>%
    #Count B in both VT and NVT
    #.[, VT := VT + B*1] %>% #count B in VT, assume it is the dominant serotype
    #.[, NVT := NVT + B*0] %>% #count B in NVT, assume it is the dominant serotype
    #.[, -"B"] %>%
    melt(measure.vars = compartments_prevalence,
         variable.name = "compartment", value.name = "modelled")
}

prepareDataForLL = function(res_prevacc, data){
  age_groups_model[to <= data[, max(to)]] %>%
    matchingAgeBreaks(data) %>% .[, age := 1:.N] %>% merge(res_prevacc, by="age") %>%
    .[, .(value = sum(value), N = sum(N)), by=c("compartment", "vaccination_group", "name.y")] %>%
    .[, .(value = sum(value), N = unique(N)), by = c("compartment", "name.y")] %>%
    .[, prev := value/N] %>%
    dcast(name.y ~ compartment, value.var = "prev") %>%
    #Count B in both VT and NVT
    #.[, VT := VT + B*1] %>% #count B in VT, assume it is the dominant serotype
    #.[, NVT := NVT + B*0] %>% #count B in NVT, assume it is the dominant serotype
    #.[, -"B"] %>%
    melt(measure.vars = compartments_prevalence,
         variable.name = "compartment", value.name = "modelled")
}

#' Calculate the log likelihood based on the pre-vaccination model output (at equilibrium)
calculateLLprevacc = function(res_prevacc, data){
  result = res_prevacc %>%
    .[grepl("idp_", population)] %>%
    aggregateModelOutput(model_params, data,
                         by_population = FALSE, by_vaccination_group = FALSE) %>%
    .[, modelled := value]
  
  #' Make sure all modelled value are positive finite proportions
  if(any(result[, modelled] < 0 | is.infinite(result[, modelled]))) return(-Inf)
  
  result %>%
    merge(data[, -c("from", "to")] %>%
            melt(measure.vars = model_params$global_settings$compartments_prevalence,
                 variable.name = "compartment", value.name = "observed"),
          by.x = c("age_group", "compartment"), by.y = c("name", "compartment")) %>%
    dcast(age_group ~ compartment, value.var = c("modelled", "observed")) %>%
    .[, .SD[, dmultinom(c(observed_S, observed_VT, observed_NVT, observed_B),
                        prob = c(modelled_S, modelled_VT, modelled_NVT, modelled_B),
                        log = TRUE)], by = c("age_group")] %>% .[, sum(V1)]
}

calculateLL_neutral = function(model_params_current, carriage_data, parallel = parallel){
  result_prevacc = runModel(model_params = model_params_current,
                            initial_state = model_params_current$global_settings$initial_states,
                            steady_state = TRUE, parallel = parallel)
  
  log_ll_total = -Inf
  #' Only continue if no errors or warning, otherwise reject sample by returning -Inf
  if(result_prevacc$status == 0 & checkModelOutput(result_prevacc$value)){
    result_prevacc = result_prevacc$value
    
    #' calculate log-likelihood of model fit to data
    log_ll_total = result_prevacc %>% calculateLLprevacc_neutral(carriage_data)
  }
  
  return(list(log_ll_total = log_ll_total))
}

#' Calculate the log likelihood based on the pre-vaccination model output (at equilibrium)
calculateLLprevacc_neutral = function(res_prevacc, data){
  result = res_prevacc %>%
    .[grepl("idp_", population)] %>%
    aggregateModelOutput(model_params, data,
                         by_population = FALSE, by_vaccination_group = FALSE) %>%
    .[, modelled := value]
  
  #' Make sure all modelled value are positive finite proportions
  if(any(result[, modelled] < 0 | is.infinite(result[, modelled]))) return(-Inf)
  
  result %>%
    merge(data[, -c("from", "to")] %>%
            melt(measure.vars = model_params$global_settings$compartments_prevalence,
                 variable.name = "compartment", value.name = "observed"),
          by.x = c("age_group", "compartment"), by.y = c("name", "compartment")) %>%
    dcast(age_group ~ compartment, value.var = c("modelled", "observed")) %>%
    .[, .SD[, dmultinom(c(observed_S, observed_VT, observed_VT2, observed_NVT, observed_NVT2, observed_B),
                        prob = c(modelled_S, modelled_VT, modelled_VT2, modelled_NVT, modelled_NVT2, modelled_B),
                        log = TRUE)], by = c("age_group")] %>% .[, sum(V1)]
}

#' Function that will be used in MCMC algorithm
calculateLL2 = function(model_params, data_unvacll, parallel = TRUE){
  res_prevacc = runModel(initial_state = model_params$state_prop_cluster, steady_state = TRUE,
                         model_params = model_params$params_unvac, parallel = parallel)
  
  #' Only continue if no errors or warning, otherwise reject sample by returning -Inf
  if(res_prevacc$status == 0){
    res_prevacc = res_prevacc$value
    if(res_prevacc[, any(value < 0) | any(value > 1)]){
      warning("Some values negative or greater than 1. Returning early.")
      return(list(log_ll_total = -Inf))
    }
    
    #' Get population size by age in all populations to calculate total prevalence for weighted mean
    N = model_params$params_unvac$trial_arms %>% seq_along %>%
      lapply(function(x, clusters) data.table(N = clusters[[x]]$parameters$N, population = names(clusters)[x]) %>% .[, age := .I],
             model_params$params_unvac$trial_arms) %>% rbindlist
    res_prevacc = res_prevacc %>% merge(N, by = c("population", "age")) %>% .[, value := value * N]
    
    #' Calculate log likelihood only for the idp population
    log_ll_total = res_prevacc[grepl("idp_", population)] %>% calculateLLprevacc(data_unvacll)
  } else {
    log_ll_total = -Inf
  }
  
  return(list(log_ll_total = log_ll_total))
}

lapplyNamed = function(x, fun, ...){
  if(is.null(names(x))) stop("Element has no names")
  #if(length(formals(fun)) != 1) stop("Function does not only have a single argument")
  if(is.list(x)){
    lapply(names(x), function(x, vals = x){ x = setNames(list(vals[[x]]), x); fun(x, ...) }, x)
  } else {
    lapply(names(x), function(x, vals = x){ x = vals[x]; fun(x, ...) }, x) 
  }
}

#' Increases or decreases the brightness of a color by a percentage of the current brightness.
#' adjust_by: A number between -1 and 1. E.g. 0.3 = 30% lighter; -0.4 = 40% darker.
#' See https://stackoverflow.com/a/54393956
adjustColBrightness = function(hex_code, adjust_by=c(-1, 0, 1)){
  if(any(adjust_by < -1) | any(adjust_by > 1)) stop("adjust_by needs to be between -1 and 1")
  
  if (nchar(hex_code) == 4) {
    hex_code = gsub("#", "", hex_code) %>%
      strsplit("") %>% .[[1]] %>% rep(each=2) %>% paste0("#", .)
  }
  
  rgb = col2rgb(hex_code)
  
  col_matrix = matrix(rep(rgb, length(adjust_by)), 3)
  adjustable_limit = col_matrix
  adjustable_limit[, which(adjust_by > 0)] = (255 - as.matrix(adjustable_limit[, which(adjust_by > 0)]))
  adjustable_limit = ceiling(adjustable_limit * matrix(rep(adjust_by, each=3), 3))
  col_matrix = col_matrix + adjustable_limit
  
  apply(col_matrix, MARGIN = 2, FUN = function(x) rgb(x[1]/255, x[2]/255, x[3]/255, 1))
}

sampleCaseCarrierRatio = function(age_groups_model, case_carrier_data){
  i = sample(case_carrier_data[, iter], 1)
  
  case_carrier_ratio_model = age_groups_model %>%
    combineAgeBreaks(case_carrier_data[iter == i] %>% dcast(from+to+name ~ st),
                     value.var = c("NVT", "VT"))
  
  return(case_carrier_ratio_model)
}

adjustCaseCarrierRatio = function(ccr, kilifi_maln_prev = 0.2, maln_RR = 2){
  base_risk = ccr/(1 - kilifi_maln_prev + kilifi_maln_prev*maln_RR)
  return(base_risk)
}

calculateLL = function(model_params_current, carriage_data, parallel = parallel){
  result_prevacc = runModel(model_params = model_params_current,
                            initial_state = model_params_current$global_settings$initial_states,
                            steady_state = TRUE, parallel = parallel)
  
  log_ll_total = -Inf
  #' Only continue if no errors or warning, otherwise reject sample by returning -Inf
  if(result_prevacc$status == 0 & checkModelOutput(result_prevacc$value)){
    result_prevacc = result_prevacc$value
    
    #' calculate log-likelihood of model fit to data
    log_ll_total = result_prevacc %>% calculateLLprevacc(carriage_data)
  }
  
  return(list(log_ll_total = log_ll_total))
}

lshtm_colours = list(black = "#000000", white = "#FFFFFF", lightgrey = "#DEDEDE", grey="#555555", darkgrey = "#AAAAAA",
               red = "#FE5000", yellow = "#FFB81C", blue = "#00AEC7", green = "#0D5257", lightgreen = "#00BF6F",
               purple = "#621244", pink = "#FFABBA")

SLIDE_WIDTH = 13.33333
SLIDE_HEIGHT = 7.5
BANNER_HEIGHT = 1.19
SLIDE_MARGIN_HORIZONTAL = 0.67
SLIDE_MARGIN_VERTICAL = 0.31
PLOT_WIDTH_ONE = (SLIDE_WIDTH - 2*SLIDE_MARGIN_HORIZONTAL)
PLOT_WIDTH_ONETWO = (SLIDE_WIDTH - 3*SLIDE_MARGIN_HORIZONTAL)/2 * 1
PLOT_WIDTH_ONETHIRD = (SLIDE_WIDTH - 3*SLIDE_MARGIN_HORIZONTAL)/3 * 1
PLOT_WIDTH_TWOTHIRD = (SLIDE_WIDTH - 3*SLIDE_MARGIN_HORIZONTAL)/3 * 2
PLOT_HEIGHT_ONE = SLIDE_HEIGHT - BANNER_HEIGHT - 2*SLIDE_MARGIN_VERTICAL
PLOT_HEIGHT_ONETWO = (SLIDE_HEIGHT - BANNER_HEIGHT - 3*SLIDE_MARGIN_VERTICAL)/2
PLOT_HEIGHT_ONETHIRD = (SLIDE_HEIGHT - BANNER_HEIGHT - 3*SLIDE_MARGIN_VERTICAL)/3 * 1
PLOT_HEIGHT_TWOTHIRD = (SLIDE_HEIGHT - BANNER_HEIGHT - 3*SLIDE_MARGIN_VERTICAL)/3 * 2
PLOT_HEIGHT_THREEFOURTH = (SLIDE_HEIGHT - BANNER_HEIGHT - 3*SLIDE_MARGIN_VERTICAL)/4*3