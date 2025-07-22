#' load/update all libraries
pacman::p_load(data.table, Rcpp, RcppArmadillo, inline, deSolve, rootSolve, magrittr, units, ggplot2, BayesianTools, ggplot2, ggh4x, lvplot, binom, tictoc)

# Can setwd by clicking on Session > Set Working Directory > To Source File Location in the toolbar in RStudio

#' set arguments if R is running interactively, otherwise arguments should be passed to the Rscript command
#'  1: working directory of pcvmr (make sure to setwd manually if running interactively)
#'  2: working directory of pcvm
#'  3: name of output folder to store objects (e.g. date, or specific scenario name)
.args = if(interactive()) c(getwd(), "../pcvm", "digaale", "simulations", 9) else commandArgs(trailingOnly = TRUE)
.args = setNames(.args, c("wd", "metavax_dir", "output_subdir", "output_simdir", "cores"))
setwd(.args["wd"])

#' Load metavax model and setup initial model_params skeleton
METAVAX_FOLDER = .args["metavax_dir"]
source(sprintf("%s/model/metavax_streppneumo_diamond.R", METAVAX_FOLDER))

#' Create folder to store output
OUTPUT_FOLDER = setOutputFolder(.args["wd"], .args["output_subdir"])
OUTPUT_SUBFOLDER = setOutputFolder(.args["wd"], sprintf("%s/%s", .args["output_subdir"], .args["output_simdir"]))

#' Set up project specific functions and parameters
source("./functions.R")
source("./data_load_all.R")
source("./model_setup.R")

# combine all out files in OUTPUT_FOLDER
out = combineOutFiles(OUTPUT_FOLDER)
set.seed(123)
posterior = samplePosterior(out, start = 6000, thin = 10, variable_names = as.character(priors$variable))

#' we need to create a function to run a single iteration of the model using a sample from the posterior,
#'  and return the correct data
#' create model params without vaccination
if(.args["output_subdir"] == "bambari"){
  model_params_prevac = modelSetup(contact_data = list(data$demography[[.args["output_subdir"]]]$displaced$contact_data$data),
                                   contacts_extra_host = data$demography[[.args["output_subdir"]]]$displaced$contacts_extra_host$data$value,
                                   contacts_prop_extra = data$demography[[.args["output_subdir"]]]$displaced$contacts_prop_extra$data$value,
                                   population_data = list(data$demography[[.args["output_subdir"]]]$displaced$population_data$data),
                                   prop_malnourished = c(data$demography[[.args["output_subdir"]]]$displaced$malnourished$data$value),
                                   migration_rate = data$demography[[.args["output_subdir"]]]$displaced$migration_rate$data$value,
                                   vaccination_groups = "unvaccinated")
  
  model_params = modelSetup(contact_data = list(data$demography[[.args["output_subdir"]]]$displaced$contact_data$data),
                            contacts_extra_host = data$demography[[.args["output_subdir"]]]$displaced$contacts_extra_host$data$value,
                            contacts_prop_extra = data$demography[[.args["output_subdir"]]]$displaced$contacts_prop_extra$data$value,
                            population_data = list(data$demography[[.args["output_subdir"]]]$displaced$population_data$data),
                            prop_malnourished = c(data$demography[[.args["output_subdir"]]]$displaced$malnourished$data$value),
                            migration_rate = data$demography[[.args["output_subdir"]]]$displaced$migration_rate$data$value)
} else {
  model_params_prevac = modelSetup(contact_data = list(data$demography[[.args["output_subdir"]]]$displaced$contact_data$data,
                                                       data$demography[[.args["output_subdir"]]]$host$contact_data$data),
                                   contacts_extra_host = data$demography[[.args["output_subdir"]]]$displaced$contacts_extra_host$data$value,
                                   contacts_prop_extra = data$demography[[.args["output_subdir"]]]$displaced$contacts_prop_extra$data$value,
                                   population_data = list(data$demography[[.args["output_subdir"]]]$displaced$population_data$data,
                                                          data$demography[[.args["output_subdir"]]]$host$population_data$data),
                                   prop_malnourished = c(data$demography[[.args["output_subdir"]]]$displaced$malnourished$data$value,
                                                         data$demography[[.args["output_subdir"]]]$host$malnourished$data$value),
                                   migration_rate = data$demography[[.args["output_subdir"]]]$displaced$migration_rate$data$value,
                                   vaccination_groups = "unvaccinated")
  
  model_params = modelSetup(contact_data = list(data$demography[[.args["output_subdir"]]]$displaced$contact_data$data,
                                                data$demography[[.args["output_subdir"]]]$host$contact_data$data),
                            contacts_extra_host = data$demography[[.args["output_subdir"]]]$displaced$contacts_extra_host$data$value,
                            contacts_prop_extra = data$demography[[.args["output_subdir"]]]$displaced$contacts_prop_extra$data$value,
                            population_data = list(data$demography[[.args["output_subdir"]]]$displaced$population_data$data,
                                                   data$demography[[.args["output_subdir"]]]$host$population_data$data),
                            prop_malnourished = c(data$demography[[.args["output_subdir"]]]$displaced$malnourished$data$value,
                                                  data$demography[[.args["output_subdir"]]]$host$malnourished$data$value),
                            migration_rate = data$demography[[.args["output_subdir"]]]$displaced$migration_rate$data$value)
}

singleRun = function(params, parallel = TRUE){
  #' parameter values used in this run
  params = c(params, additional_posteriors$sampler())
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
      model_params_current = updateParameters(params, model_params)
    } else {
      model_params_current = model_params_current %>%
        setParameter(level = "vaccination_group", population = c("idp_malnourished", "idp_non_malnourished"),
                     vaccination_group = "unvaccinated",
                     key = "coverage_c",
                     value = list(
                       list(value = getVaccineCoverage(model_params_current$global_settings$age_groups_model,
                                                       c(set_units(6, "weeks"), set_units(upper_age_limit, "years")), 0.85),
                            time = 0, coverage_to = c(rep("pcv1", age_groups_model[to <= set_units(12, "months"), .N]),
                                                      rep("pcv2", age_groups_model[from >= set_units(12, "months"), .N]))),
                       #all infants receive a 2nd dose and move to pcv2
                       list(value = getVaccineCoverage(model_params_current$global_settings$age_groups_model,
                                                       c(set_units(6, "weeks"), set_units(1, "years")), 1.00),
                            time = 1, coverage_to = c(rep("pcv2", age_groups_model[to <= set_units(12, "months"), .N]),
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

#' iterations to run
iterations = 500
parameter_run_values = posterior[sample(seq_len(.N), iterations), -"chain"] %>%
  cbind(additional_posteriors$sampler(iterations))

tictoc::tic()
posterior_runs = parallelModelRuns(singleRun,
                                   parameter_run_values,
                                   cores = as.numeric(.args["cores"]))
tictoc::toc()

saveRDS(posterior_runs, sprintf("%s/posterior_runs_2dose_infants.RDS", OUTPUT_SUBFOLDER))
saveRDS(parameter_run_values, sprintf("%s/parameter_values_2dose_infants.RDS", OUTPUT_SUBFOLDER))
