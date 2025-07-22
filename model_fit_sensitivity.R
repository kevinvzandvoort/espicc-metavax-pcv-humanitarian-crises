pacman::p_load(data.table, Rcpp, RcppArmadillo, inline, deSolve, rootSolve, readxl, magrittr, binom, qs, units, BayesianTools)

.args = if(interactive()) c(getwd(), "../pcvm", "contacts_host", 0, 1, 1) else commandArgs(trailingOnly = TRUE)
.args = setNames(.args, c("wd", "metavax_dir", "sens_var", "sens_val", "chain", "i"))
setwd(.args["wd"])

#' Set seed for RNG
set.seed(as.numeric(.args["chain"]))

#' Load metavax model and setup initial model_params skeleton
METAVAX_FOLDER = .args["metavax_dir"]
source(sprintf("%s/model/metavax_streppneumo_diamond.R", METAVAX_FOLDER))

if(!.args["sens_var"] %in% c("contacts_host", "migration_rate", "mal_trans"))
  stop("variable not implemented")

#' Create folder to store output
if(!dir.exists("output/sensitivity")) dir.create("output/sensitivity")
if(!dir.exists(sprintf("output/sensitivity/%s", .args["sens_var"]))) dir.create(sprintf("output/sensitivity/%s", .args["sens_var"]))
OUTPUT_FOLDER = setOutputFolder(.args["wd"], sprintf("sensitivity/%s/%s", .args["sens_var"], .args["sens_val"]))

#' Set up project specific functions and parameters
source("./functions.R")
source("./data_load_all.R")
source("./model_setup.R")

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

#' create model params without vaccination
model_params = modelSetup(contacts_extra_host = sens_contacts_host,
                          migration_rate = sens_migration_rate,
                          malnourished_rr_transmission = sens_maln_trans,
                          vaccination_groups = "unvaccinated")
#' Create log-likelihood function for BayesianTools
ll = function(params, parallel = TRUE){
  names(params) = priors[, variable]
  params = c(params, additional_posteriors$sampler())
  
  #make sure updateParameters reads named VECTOR
  model_params_current = updateParameters(params, model_params)
  out = calculateLL(model_params_current, data$epidemiology$prevalence$digaale_pneumosil$data, parallel = parallel)
  
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
