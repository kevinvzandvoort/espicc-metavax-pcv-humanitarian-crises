#' load/update all libraries
pacman::p_load(data.table, Rcpp, RcppArmadillo, inline, deSolve, rootSolve, magrittr, units, ggplot2, BayesianTools, ggplot2, ggh4x, lvplot, binom)

# Can setwd by clicking on Session > Set Working Directory > To Source File Location in the toolbar in RStudio

#' set arguments if R is running interactively, otherwise arguments should be passed to the Rscript command
#'  1: working directory of pcvmr (make sure to setwd manually if running interactively)
#'  2: working directory of pcvm
#'  3: name of output folder to store objects (e.g. date, or specific scenario name)
.args = if(interactive()) c(getwd(), "../pcvm", "digaale") else commandArgs(trailingOnly = TRUE)
.args = setNames(.args, c("wd", "metavax_dir", "output_subdir"))
setwd(.args["wd"])

#' Load metavax model and setup initial model_params skeleton
METAVAX_FOLDER = .args["metavax_dir"]
if(grepl("neutral", .args["output_subdir"])){
  source(sprintf("%s/model/metavax_streppneumo_neutral.R", METAVAX_FOLDER))
} else {
  source(sprintf("%s/model/metavax_streppneumo_diamond.R", METAVAX_FOLDER))  
}

#' Create folder to store output
OUTPUT_FOLDER = setOutputFolder(.args["wd"], .args["output_subdir"])

#' Set up project specific functions and parameters
source("./functions.R")
source("./data_load_all.R")
source("./model_setup.R")

# combine all out files in OUTPUT_FOLDER
out = combineOutFiles(OUTPUT_FOLDER)

#' check trace plot
#plot(out)
burnin = 6000
#plot(out, start=burnin)

#' create geglman Diagnostics
#gelmanDiagnostics(out, plot=TRUE, start = burnin)
#correlationPlot(out, start = burnin)

#' Set seed for RNG
set.seed(123)

gd = gelmanDiagnostics(out, start = burnin)
gelman_table = data.table(variable = c(priors$variable, "Multivariate PSRF"),
                          r_hat_est = round(c(gd$psrf[, 1], gd$mpsrf), 3),
                          r_hat_upper = round(c(gd$psrf[, 2], NA), 3))
fwrite(gelman_table, sprintf("%s/gelman_table.csv", OUTPUT_FOLDER))

#' sample from the posterior
posterior = samplePosterior(out, start = burnin, thin = 10, variable_names = as.character(priors$variable))

coda::crosscorr(coda::mcmc(posterior[, -"chain"])) %>% round(2) %>%
  fwrite(sprintf("%s/cross_correlation.csv", OUTPUT_FOLDER))

#' Let's plot the trace plots in a more readable form
createTracePlot(posterior)
ggsave(sprintf("%s/posterior_trace.png", OUTPUT_FOLDER), width=10, height=8)

#' Let's also compare the posterior distributions to the prior distributions
createPriorPosteriorPlot(priors, posterior)
ggsave(sprintf("%s/posterior_density.png", OUTPUT_FOLDER), width=12, height=8)

#' we need to create a function to run a single iteration of the model using a sample from the posterior,
#'  and return the correct data
#' create model params without vaccination
setting = strsplit(.args["output_subdir"], "_", TRUE)[[1]][1]
if(setting == "bambari"){
  model_params = modelSetup(contact_data = list(data$demography[[setting]]$displaced$contact_data$data),
                            contacts_extra_host = data$demography[[setting]]$displaced$contacts_extra_host$data$value,
                            contacts_prop_extra = data$demography[[setting]]$displaced$contacts_prop_extra$data$value,
                            population_data = list(data$demography[[setting]]$displaced$population_data$data),
                            prop_malnourished = c(data$demography[[setting]]$displaced$malnourished$data$value),
                            migration_rate = data$demography[[setting]]$displaced$migration_rate$data$value,
                            vaccination_groups = "unvaccinated")
} else {
  model_params = modelSetup(contact_data = list(data$demography[[setting]]$displaced$contact_data$data,
                                                data$demography[[setting]]$host$contact_data$data),
                            contacts_extra_host = data$demography[[setting]]$displaced$contacts_extra_host$data$value,
                            contacts_prop_extra = data$demography[[setting]]$displaced$contacts_prop_extra$data$value,
                            population_data = list(data$demography[[setting]]$displaced$population_data$data,
                                                   data$demography[[setting]]$host$population_data$data),
                            prop_malnourished = c(data$demography[[setting]]$displaced$malnourished$data$value,
                                                  data$demography[[setting]]$host$malnourished$data$value),
                            migration_rate = data$demography[[setting]]$displaced$migration_rate$data$value,
                            vaccination_groups = "unvaccinated")
}

#' set initial states for neutral model
if(grepl("neutral", .args["output_subdir"])){
  vaccination_groups = "unvaccinated"
  model_params = model_params %>%
    setParameter(level = "global", key = "initial_states",
                 value = c(0.8, 0.1, 0.1, 0, 0, 0) %>%
                   rep(each = model_params$global_settings$age_groups_model[, .N]) %>%
                   c(c(0, 0, 0, 0, 0, 0) %>% rep(each = model_params$global_settings$age_groups_model[, .N]) %>%
                       rep(length(vaccination_groups) - 1)) %>%
                   rep(length(model_params$populations)))
}

singleRun = function(params, parallel = TRUE){
  model_params_current = updateParameters(params, model_params)
  result_prevacc = runModel(model_params = model_params_current,
                            initial_state = model_params_current$global_settings$initial_states,
                            steady_state = TRUE, parallel = parallel)
  result_prevacc = result_prevacc$value
  result_aggregated = result_prevacc %>%
    .[grepl("idp_", population)] %>%
    aggregateModelOutput(model_params,
                         data$epidemiology$prevalence$digaale_pneumosil$data,
                         by_population = FALSE, by_vaccination_group = FALSE)
  
  return(result_aggregated)
}

#' iterations to run
iterations = 500
parameter_run_values = posterior[sample(seq_len(.N), iterations), -"chain"] %>%
  cbind(additional_posteriors$sampler(iterations))

#' test
singleRun(unlist(parameter_run_values[1, ]))

tictoc::tic()
posterior_runs = parallelModelRuns(singleRun,
                                   parameter_run_values,
                                   cores = 5)
tictoc::toc()

#' set initial states for neutral model
if(grepl("neutral", .args["output_subdir"])){
  data_observed = data$epidemiology$prevalence$digaale_pneumosil_alternative$data %>%
    .[, N := S + NVT + NVT2 + VT + VT2 + B, by="name"] %>%
    .[, age_group := name] %>%
    .[, -"name"] %>%
    melt(measure.vars = model_params$global_settings$compartments_prevalence, variable.name="compartment") %>%
    .[, c("mid", "low95", "high95") := binom.confint(value, N, methods="exact")[, c("mean", "lower", "upper")]]
} else {
  data_observed = data$epidemiology$prevalence$digaale_pneumosil$data %>%
    .[, N := S + NVT + VT + B, by="name"] %>%
    .[, age_group := name] %>%
    .[, -"name"] %>%
    melt(measure.vars = model_params$global_settings$compartments_prevalence, variable.name="compartment") %>%
    .[, c("mid", "low95", "high95") := binom.confint(value, N, methods="exact")[, c("mean", "lower", "upper")]]
}

posterior_runs %>%
  ggplot(aes(x=age_group, y=value, fill = compartment, colour = compartment))+
  facet_grid(.~compartment)+
  geom_segment(data = data_observed, aes(x = age_group, xend = age_group, y=low95, yend=high95),
               colour=adjustColBrightness(lshtm_colours$grey, 0.8), size=6)+
  geom_errorbarh(data = data_observed, aes(xmin = age_group, xmax = age_group, y=mid, colour=compartment),
                 colour=adjustColBrightness(lshtm_colours$black, 0.3), linewidth=6, height=0.02)+
  geom_lv(varwidth = TRUE, width.method = "height", aes(alpha = after_stat(LV)),
          outlier.colour = adjustColBrightness(lshtm_colours$grey, 0.25), outlier.shape = 19, outlier.size = 1,
          outlier.stroke = 0.25, percent=5)+
  scale_fill_manual(values = c("S" = lshtm_colours$lightgreen,
                               "VT" = lshtm_colours$red, 
                               "NVT" =  lshtm_colours$blue,
                               "VT2" = adjustColBrightness(lshtm_colours$red, 0.4),
                               "B"  = lshtm_colours$yellow,
                               "NVT2" = adjustColBrightness(lshtm_colours$blue, 0.4)))+
  scale_colour_manual(values = c("S" = lshtm_colours$lightgreen,
                                 "VT" = lshtm_colours$red, 
                                 "NVT" =  lshtm_colours$blue,
                                 "VT2" = adjustColBrightness(lshtm_colours$red, 0.4),
                                 "B"  = lshtm_colours$yellow,
                                 "NVT2" = adjustColBrightness(lshtm_colours$blue, 0.4)))+
  scale_alpha_manual(values = c(1, seq(0.8, 0.2, length.out=5)))+
  guides(alpha = "none")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  scale_y_continuous(labels=scales::percent, expand=c(0.01, 0, 0, 0), lim=c(0, 1))+
  labs(x="Age group", y="Prevalence", colour="Compartment", fill="Compartment")
ggsave(sprintf("%s/model_fit.png", OUTPUT_FOLDER), width=10, height=3)
