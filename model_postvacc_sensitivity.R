pacman::p_load(data.table, Rcpp, RcppArmadillo, inline, deSolve, rootSolve, readxl, magrittr, binom, qs, units, BayesianTools)

.args = if(interactive()) c(getwd(), "../pcvm", "vac_cov", 0.85, "simulations", 5) else commandArgs(trailingOnly = TRUE)
.args = setNames(.args, c("wd", "metavax_dir", "sens_var", "sens_val", "output_simdir", "cores"))
setwd(.args["wd"])

#' Set seed for RNG
set.seed(123)

#' Load metavax model and setup initial model_params skeleton
METAVAX_FOLDER = .args["metavax_dir"]
source(sprintf("%s/model/metavax_streppneumo_diamond.R", METAVAX_FOLDER))

if(!.args["sens_var"] %in% c("contacts_host", "migration_rate", "mal_trans", "vac_cov", "vac_eff", "vac_eff_dur"))
  stop("variable not implemented")

#' Create folder to store output
if(!dir.exists("output/sensitivity")) dir.create("output/sensitivity")
if(!dir.exists(sprintf("output/sensitivity/%s", .args["sens_var"]))) dir.create(sprintf("output/sensitivity/%s", .args["sens_var"]))
OUTPUT_FOLDER = setOutputFolder(.args["wd"], sprintf("sensitivity/%s/%s", .args["sens_var"], .args["sens_val"]))
OUTPUT_SUBFOLDER = setOutputFolder(OUTPUT_FOLDER, .args["output_simdir"])

#' Set up project specific functions and parameters
source("./functions.R")
source("./data_load_all.R")
source("./model_setup.R")

# combine all out files
if(.args["sens_var"] %in% c("contacts_host", "migration_rate", "mal_trans")){
  out = combineOutFiles(OUTPUT_FOLDER)
  posterior = samplePosterior(out, start = 6000, thin = 10, variable_names = as.character(priors$variable))
} else {
  # use main Digaale fit for sensitivity analyses that did not require refitting
  out = combineOutFiles(sprintf("%s/%s", .args["wd"], "output/digaale"))
  posterior = samplePosterior(out, start = 6000, thin = 10, variable_names = as.character(priors$variable)) 
}

#' iterations to run
set.seed(123)
iterations = 200
parameter_run_values = posterior[sample(seq_len(.N), iterations), -"chain"] %>%
  cbind(additional_posteriors$sampler(iterations))

#' set parameter value
if(.args["sens_var"] == "contacts_host"){
  sens_contacts_host = as.numeric(.args["sens_val"])
} else {
  sens_contacts_host = data$demography$digaale$displaced$contacts_extra_host$data$value
}
if(.args["sens_var"] == "migration_rate"){
  #' parameterize it by average number of years until people migrate
  sens_migration_rate = as.numeric(.args["sens_val"])
  if(sens_migration_rate > 0) sens_migration_rate = 1/(sens_migration_rate * 365)
} else {
  sens_migration_rate = data$demography$digaale$displaced$migration_rate$data$value
}
if(.args["sens_var"] == "mal_trans"){
  sens_maln_trans = as.numeric(.args["sens_val"])
} else {
  sens_maln_trans = data$epidemiology$malnourished_RR_transmission$data$value
}
if(.args["sens_var"] == "vac_cov"){
  sens_vac_cov = as.numeric(.args["sens_val"])
} else {
  sens_vac_cov = 0.85
}
if(.args["sens_var"] == "vac_eff"){
  parameter_run_values[, efficacy_carriage := as.numeric(.args["sens_val"])]
}
if(.args["sens_var"] == "vac_eff_dur"){
  parameter_run_values[, efficacy_duration := as.numeric(.args["sens_val"])]
}

#' create model params without and with vaccination
model_params_prevac = modelSetup(contacts_extra_host = sens_contacts_host,
                                 migration_rate = sens_migration_rate,
                                 malnourished_rr_transmission = sens_maln_trans,
                                 vaccination_groups = "unvaccinated")

model_params = modelSetup(contacts_extra_host = sens_contacts_host,
                          migration_rate = sens_migration_rate,
                          malnourished_rr_transmission = sens_maln_trans)

singleRun = function(params, parallel = TRUE){
  #' parameter values used in this run
  model_params_prevac_current = updateParameters(params, model_params_prevac)
  model_params_current = updateParameters(params, model_params)
  
  #' first run until equilibrium
  result_prevacc = runModel(model_params = model_params_prevac_current,
                            initial_state = model_params_prevac_current$global_settings$initial_states,
                            steady_state = TRUE, parallel = parallel)
  result_prevacc = result_prevacc$value
  prevacc_equilibrium = result_prevacc %>% eqStatesVaccinate2(model_params_current) %>% .[, value]
  
  #' run with vaccination
  upper_age_limits = c(0, 1, 2, 5, 10, 15)
  results_postvacc = lapply(upper_age_limits, function(upper_age_limit){
    message(sprintf("Simulating <%s campaign", upper_age_limit))
    
    if(upper_age_limit == 0){
      model_params_current = updateParameters(params, model_params) %>%
        setParameter(level = "vaccination_group", population = c("idp_malnourished", "idp_non_malnourished"),
                     vaccination_group = "unvaccinated",
                     key = "coverage_c",
                     value = list(
                       list(value = getVaccineCoverage(model_params_current$global_settings$age_groups_model,
                                                       set_units(0, "months"), 0),
                            time = 0, coverage_to = c(rep("pcv1", age_groups_model[, .N]))))) %>%
        renameCoverageTo2()
    } else {
      model_params_current = model_params_current %>%
        setParameter(level = "vaccination_group", population = c("idp_malnourished", "idp_non_malnourished"),
                     vaccination_group = "unvaccinated",
                     key = "coverage_c",
                     value = list(
                       list(value = getVaccineCoverage(model_params_current$global_settings$age_groups_model,
                                                       c(set_units(6, "weeks"), set_units(upper_age_limit, "years")), sens_vac_cov),
                            time = 0, coverage_to = c(rep("pcv1", age_groups_model[to <= set_units(12, "months"), .N]),
                                                      rep("pcv2", age_groups_model[from >= set_units(12, "months"), .N]))))) %>%
        renameCoverageTo2()  
    }
    
    result_postvacc = runModel(model_params = model_params_current,
                               initial_state = prevacc_equilibrium,
                               steady_state = FALSE,
                               times = seq(0, 365 * 5, ifelse(model_params_current$global_settings$model_solver_difference,
                                                              model_params_current$global_settings$solver_difference_delta_t,
                                                              1)),
                               incidence = TRUE)
    result_postvacc = result_postvacc$value
    result_postvacc = result_postvacc %>% 
      processIncidence(model_params_current)
    result_postvacc = result_postvacc[grepl("idp_", population)]
    result_postvacc_infants = result_postvacc %>%
      aggregateModelOutput(model_params_current,
                           setAgeBreaks(1),
                           by_population = FALSE, by_vaccination_group = FALSE) %>%
      .[age_group == "[0y, 1y)"]
    result_postvacc_all = result_postvacc %>%
      aggregateModelOutput(model_params_current,
                           setAgeBreaks(0),
                           by_population = FALSE, by_vaccination_group = FALSE)
    
    result_postvacc = rbind(result_postvacc_infants,
                            result_postvacc_all)
    
    result_postvacc[, vaccine_strategy := factor(sprintf("<%s", upper_age_limit),
                                                 sprintf("<%s", upper_age_limits))]
    
    return(result_postvacc)
  })
  
  results_postvacc = rbindlist(results_postvacc)
  
  #' only save required values
  results_postvacc = results_postvacc[outcome == "prevalence" & compartment %in% c("VT", "B")] %>%
    .[, .(value = sum(value), compartment = "VT"), by = c("age_group", "outcome", "time", "vaccine_strategy")] %>%
    rbind(results_postvacc[outcome == "incidence", .(value = sum(value), compartment = "all"), by = c("age_group", "outcome", "time", "vaccine_strategy")])
  
  return(results_postvacc)
}

tictoc::tic()
posterior_runs = parallelModelRuns(singleRun,
                                   parameter_run_values,
                                   cores = as.numeric(.args["cores"]))
tictoc::toc()

saveRDS(posterior_runs, sprintf("%s/posterior_runs_main_v2.RDS", OUTPUT_SUBFOLDER))
saveRDS(parameter_run_values, sprintf("%s/parameter_run_values_main_v2.RDS", OUTPUT_SUBFOLDER))

