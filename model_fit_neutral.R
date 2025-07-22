pacman::p_load(data.table, Rcpp, RcppArmadillo, inline, deSolve, rootSolve, readxl, magrittr, binom, qs, units, BayesianTools)

.args = if(interactive()) c(getwd(), "../pcvm", "digaale", 2, 1) else commandArgs(trailingOnly = TRUE)
.args = setNames(.args, c("wd", "metavax_dir", "output_subdir", "chain", "i"))
setwd(.args["wd"])

#' Set seed for RNG
set.seed(as.numeric(.args["chain"]))

#' Load metavax model and setup initial model_params skeleton
METAVAX_FOLDER = .args["metavax_dir"]
source(sprintf("%s/model/metavax_streppneumo_neutral.R", METAVAX_FOLDER))

#' Create folder to store output
OUTPUT_FOLDER = setOutputFolder(.args["wd"], .args["output_subdir"])

#' Set up project specific functions and parameters
source("./functions.R")
source("./data_load_all.R")
source("./model_setup.R")

#' create model params without vaccination
if(.args["output_subdir"] == "bambari"){
  model_params = modelSetup(contact_data = list(data$demography[[.args["output_subdir"]]]$displaced$contact_data$data),
                            contacts_extra_host = data$demography[[.args["output_subdir"]]]$displaced$contacts_extra_host$data$value,
                            contacts_prop_extra = data$demography[[.args["output_subdir"]]]$displaced$contacts_prop_extra$data$value,
                            population_data = list(data$demography[[.args["output_subdir"]]]$displaced$population_data$data),
                            prop_malnourished = c(data$demography[[.args["output_subdir"]]]$displaced$malnourished$data$value),
                            migration_rate = data$demography[[.args["output_subdir"]]]$displaced$migration_rate$data$value,
                            vaccination_groups = "unvaccinated")
} else {
  model_params = modelSetup(contact_data = list(data$demography[[.args["output_subdir"]]]$displaced$contact_data$data,
                                                data$demography[[.args["output_subdir"]]]$host$contact_data$data),
                            contacts_extra_host = data$demography[[.args["output_subdir"]]]$displaced$contacts_extra_host$data$value,
                            contacts_prop_extra = data$demography[[.args["output_subdir"]]]$displaced$contacts_prop_extra$data$value,
                            population_data = list(data$demography[[.args["output_subdir"]]]$displaced$population_data$data,
                                                   data$demography[[.args["output_subdir"]]]$host$population_data$data),
                            prop_malnourished = c(data$demography[[.args["output_subdir"]]]$displaced$malnourished$data$value,
                                                  data$demography[[.args["output_subdir"]]]$host$malnourished$data$value),
                            migration_rate = data$demography[[.args["output_subdir"]]]$displaced$migration_rate$data$value,
                            vaccination_groups = "unvaccinated")
}

#' set initial states for neutral model
vaccination_groups = "unvaccinated"
model_params = model_params %>%
  setParameter(level = "global", key = "initial_states",
             value = c(0.8, 0.1, 0.1, 0, 0, 0) %>%
               rep(each = model_params$global_settings$age_groups_model[, .N]) %>%
               c(c(0, 0, 0, 0, 0, 0) %>% rep(each = model_params$global_settings$age_groups_model[, .N]) %>%
                   rep(length(vaccination_groups) - 1)) %>%
               rep(length(model_params$populations)))

#' Create log-likelihood function for BayesianTools
ll = function(params, parallel = TRUE){
  names(params) = priors[, variable]
  params = c(params, additional_posteriors$sampler())
  
  #make sure updateParameters reads named VECTOR
  model_params_current = updateParameters(params, model_params)
  out = calculateLL_neutral(model_params_current, data$epidemiology$prevalence$digaale_pneumosil_alternative$data, parallel = parallel)
  
  return(out$log_ll_total)
}

#' Test ll function
if(interactive()){
  testLL(ll, prior)
}

#' fit model using BayesianTools
bayesianSetup = createBayesianSetup(likelihood = ll, prior = prior, parallel = 3)
settings = list(burnin = 100*3, iterations = 1000*3, nrChains = 1, consoleUpdates = 10)
fitBT(bayesianSetup, settings, OUTPUT_FOLDER, as.numeric(.args["chain"]), as.numeric(.args["i"]))
