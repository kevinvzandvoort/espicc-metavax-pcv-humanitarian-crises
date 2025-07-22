#' load/update all libraries
pacman::p_load(data.table, Rcpp, RcppArmadillo, inline, deSolve, rootSolve, magrittr, units, ggplot2, BayesianTools, ggplot2, ggh4x, lvplot, binom, tictoc, patchwork)

# Can setwd by clicking on Session > Set Working Directory > To Source File Location in the toolbar in RStudio

#' set arguments if R is running interactively, otherwise arguments should be passed to the Rscript command
#'  1: working directory of pcvmr (make sure to setwd manually if running interactively)
#'  2: working directory of pcvm
#'  3: name of output folder to store objects (e.g. date, or specific scenario name)
.args = if(interactive()) c(getwd(), "../pcvm", "digaale", "simulations") else commandArgs(trailingOnly = TRUE)
.args = setNames(.args, c("wd", "metavax_dir", "output_subdir", "output_simdir"))
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
OUTPUT_SUBFOLDER = setOutputFolder(.args["wd"], sprintf("%s/%s", .args["output_subdir"], .args["output_simdir"]))
OUTPUT_RESULTSFOLDER = setOutputFolder(.args["wd"], sprintf("%s/%s/%s", .args["output_subdir"], .args["output_simdir"], "results"))

#' Set up project specific functions and parameters
source("./functions.R")
source("./data_load_all.R")
source("./model_setup.R")

#' create model params without vaccination
setting = strsplit(.args["output_subdir"], "_", TRUE)[[1]][1]
if(setting == "maiduguri") setting = "maiduguri_acute"
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

if(grepl("neutral", .args["output_subdir"])){
  #' set initial states for neutral model
  vaccination_groups = "unvaccinated"
  model_params_prevac = model_params_prevac %>%
    setParameter(level = "global", key = "initial_states",
                 value = c(0.8, 0.1, 0.1, 0, 0, 0) %>%
                   rep(each = model_params$global_settings$age_groups_model[, .N]) %>%
                   c(c(0, 0, 0, 0, 0, 0) %>% rep(each = model_params$global_settings$age_groups_model[, .N]) %>%
                       rep(length(vaccination_groups) - 1)) %>%
                   rep(length(model_params$populations)))
}

vaccination_strategies_plot_settings = list(
  list(model_name = "<0", name = "No vaccination", target = T, dosing = F, colour = lshtm_colours$black),
  list(model_name = "<1", name = "<1y", target = T, dosing = T, colour = lshtm_colours$red),
  list(model_name = "<2", name = "<2y", target = T, dosing = T, colour = lshtm_colours$yellow),
  list(model_name = "<5", name = "<5y", target = T, dosing = T, colour = lshtm_colours$blue),
  list(model_name = "<10", name = "<10y", target = T, dosing = T, colour = lshtm_colours$green),
  list(model_name = "<15", name = "<15y", target = T, dosing = T, colour = lshtm_colours$lightgreen)) %>%
  lapply(as.data.table) %>% rbindlist()
vaccination_strategies_plot_settings[, name := factor(name, vaccination_strategies_plot_settings$name)]

posterior_runs_noroutine = readRDS(sprintf("%s/posterior_runs_main_v2.RDS", OUTPUT_SUBFOLDER))
posterior_runs_noroutine[, routine_coverage := 0]

posterior_runs = readRDS(sprintf("%s/posterior_runs_routine.RDS", OUTPUT_SUBFOLDER))
posterior_runs_parameter_values = readRDS(sprintf("%s/parameter_run_values_routine.RDS", OUTPUT_SUBFOLDER))
posterior_runs_parameter_values[, run := .I]

posterior_runs = posterior_runs %>% merge(posterior_runs_parameter_values[, c("run", "routine_coverage")], by = "run")
posterior_runs = posterior_runs %>% rbind(posterior_runs_noroutine)

routine_coverage_labels = c("0" = "0%", "0.2" = "20%", "0.4" = "40%", "0.6" = "60%", "0.8" = "80%")

posterior_runs_summarized = posterior_runs %>%
  .[, .(med = median(value),
        low50=quantile(value, 0.25), high50=quantile(value, 0.75),
        low95=quantile(value, 0.025), high95=quantile(value, 0.975)),
    by=c("age_group", "vaccine_strategy", "outcome", "time", "routine_coverage")]

(plot_prevalence_fig = posterior_runs_summarized[outcome == "prevalence"] %>%
  ggplot(aes(x=time, y=med,
             colour=factor(vaccine_strategy, vaccination_strategies_plot_settings$model_name, vaccination_strategies_plot_settings$name),
             fill=factor(vaccine_strategy, vaccination_strategies_plot_settings$model_name, vaccination_strategies_plot_settings$name),
             group = paste0(vaccine_strategy, age_group)))+
  facet_grid(factor(age_group, c("[0y, 120y)", "[0y, 1y)"), c("all ages", "infants"))~factor(routine_coverage, names(routine_coverage_labels), routine_coverage_labels), scales="free")+
  geom_ribbon(alpha=0.2, aes(ymin=low95, ymax=high95, colour = NULL))+
  geom_line(linewidth=2)+
  scale_colour_manual(values = vaccination_strategies_plot_settings[, c("name", "colour")] %>% as.matrix("name") %>% .[,1])+
  scale_fill_manual(values = vaccination_strategies_plot_settings[, c("name", "colour")] %>% as.matrix("name") %>% .[,1])+
  theme_minimal()+
  labs(x="Years since PCV campaign", y="VT carriage prevalence", colour="Vaccine strategy", fill="Vaccine strategy")+
  theme(legend.position="bottom", panel.grid.minor.y = element_blank(), strip.text = element_text(face="bold", size=14),
        axis.title = element_text(size=14), axis.text = element_text(size = 12), legend.title = element_text(size = 14),
        legend.text = element_text(size = 12), panel.spacing.y=unit(5, "mm"))+
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5) * 365, labels = paste0(c("0", 1, 2, 3, 4, 5), "y"))+
  scale_y_continuous(labels = scales::percent, lim=c(0, NA))+
  geom_vline(data=data.table(age_group = factor("[0y, 1y)", c("[0y, 120y)", "[0y, 1y)")), time = 1*365),
             aes(xintercept=time, x=NULL, colour=NULL, fill=NULL, group=NULL),
             colour="#000000", linewidth=1, linetype=2))

ggsave(plot_prevalence_fig, filename = sprintf("%s/figure_prevalence_routine.png", OUTPUT_RESULTSFOLDER),
       width = PLOT_WIDTH_ONE, height = PLOT_HEIGHT_ONE, units = "in", bg = "#FFFFFF")

posterior_runs_summarized_impact = posterior_runs[outcome == "incidence"] %>%
  dcast(...~vaccine_strategy, value.var = "value") %>%
  melt(measure.vars = vaccination_strategies_plot_settings[dosing == TRUE, model_name],
       variable.name = "vaccine_strategy") %>%
  .[, impact := 1 - value/`<0`] %>%
  .[, .(med = median(impact),
        low50=quantile(impact, 0.25), high50=quantile(impact, 0.75),
        low95=quantile(impact, 0.025), high95=quantile(impact, 0.975)),
    by=c("age_group", "vaccine_strategy", "outcome", "time", "routine_coverage")]

(plot_incidence_impact = posterior_runs_summarized_impact %>%
  ggplot(aes(x=time, y=med,
             colour=factor(vaccine_strategy, vaccination_strategies_plot_settings$model_name, vaccination_strategies_plot_settings$name),
             fill=factor(vaccine_strategy, vaccination_strategies_plot_settings$model_name, vaccination_strategies_plot_settings$name),
             group = paste0(vaccine_strategy, age_group)))+
    facet_grid(factor(age_group, c("[0y, 120y)", "[0y, 1y)"), c("all ages", "infants"))~factor(routine_coverage, names(routine_coverage_labels), routine_coverage_labels), scales="free")+
  geom_ribbon(alpha=0.2, aes(ymin=low95, ymax=high95, colour = NULL))+
  geom_line(size=2)+
  scale_colour_manual(values = vaccination_strategies_plot_settings[, c("name", "colour")] %>% as.matrix("name") %>% .[,1])+
  scale_fill_manual(values = vaccination_strategies_plot_settings[, c("name", "colour")] %>% as.matrix("name") %>% .[,1])+
  theme_minimal()+labs(x="Years since PCV campaign", y="Daily impact on severe\npneumococcal disease cases", colour="Vaccine strategy", fill="Vaccine strategy")+
  theme(legend.position="bottom", panel.grid.minor.y = element_blank(), strip.text = element_text(face="bold", size=14),
        axis.title = element_text(size=14), axis.text = element_text(size = 12), legend.title = element_text(size = 14),
        legend.text = element_text(size = 12), panel.spacing.y=unit(5, "mm"))+
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5) * 365, labels = paste0(c("0", 1, 2, 3, 4, 5), "y"))+
  scale_y_continuous(labels = scales::percent, lim=c(0, NA))+
  geom_vline(data=data.table(age_group = factor("[0y, 1y)", c("[0y, 120y)", "[0y, 1y)")),
                             routine_coverage = as.numeric(names(routine_coverage_labels)),
                             time = 1*365),
             aes(xintercept=time, x=NULL, colour=NULL, fill=NULL, group=NULL),
             colour="#000000", linewidth=1, linetype=2))

ggsave(plot_incidence_impact, filename = sprintf("%s/figure_incidence_routine.png", OUTPUT_RESULTSFOLDER),
       width = PLOT_WIDTH_ONE, height = PLOT_HEIGHT_ONE, units = "in", bg = "#FFFFFF")

population_size_all = getParameter(model_params, key = "N", level = "population", population = c("idp_non_malnourished", "idp_malnourished")) %>%
  merge(model_params$global_settings$age_groups_model, by = "age") %>%
  .[, .(N = sum(N)), by = c("age", "name", "from", "to")] %>%
  combineAgeBreaks(x = setAgeBreaks(0), y = ., method = "sum", value.var = "N")
  
population_size_infants = getParameter(model_params, key = "N", level = "population", population = c("idp_non_malnourished", "idp_malnourished")) %>%
  merge(model_params$global_settings$age_groups_model, by = "age") %>%
  .[, .(N = sum(N)), by = c("age", "name", "from", "to")] %>%
  combineAgeBreaks(x = setAgeBreaks(1), y = ., method = "sum", value.var = "N") %>%
  .[to <= set_units(1, "year")]

population_size = rbind(population_size_all, population_size_infants)

upper_age_limits = c(1, 2, 5, 10, 15)
vaccine_doses = lapply(upper_age_limits, function(upper_age_limit){
  getParameter(model_params, key = "N", level = "population", population = c("idp_non_malnourished", "idp_malnourished")) %>%
    merge(model_params$global_settings$age_groups_model, by = "age") %>%
    .[, .(N = sum(N)), by = c("age", "name", "from", "to")] %>%
    .[, coverage := getVaccineCoverage(model_params$global_settings$age_groups_model,
                                       c(set_units(6, "weeks"), set_units(upper_age_limit, "year")), 0.85)] %>%
    .[, .(vaccine_strategy = sprintf("<%s", upper_age_limit), doses = sum(N * coverage))]
}) %>% rbindlist()

upper_time_limits = c(0.5, 1, 2, 3)
incidence_impact = lapply(upper_time_limits, function(upper_time_limit){
  posterior_runs %>%
    .[outcome == "incidence"] %>%
    .[time <= 365 * upper_time_limit] %>%
    .[, .(value = sum(value)), by = c("age_group", "vaccine_strategy", "run", "routine_coverage")] %>%
    dcast(...~vaccine_strategy, value.var = "value") %>%
    melt(measure.vars = vaccination_strategies_plot_settings[dosing == TRUE, model_name],
         variable.name = "vaccine_strategy", value.name = "value") %>%
    .[, .(absolute_impact = (`<0` - value) * 5.98,
          relative_impact = 1 - value/`<0`), by = c("age_group", "vaccine_strategy", "run", "routine_coverage")] %>%
    merge(population_size[, c("name", "N")], by.y = "name", by.x = "age_group") %>%
    merge(vaccine_doses, by = "vaccine_strategy") %>%
    .[, NNV := doses/absolute_impact] %>%
    .[, absolute_impact_p10000 := absolute_impact/N * 10000] %>%
    .[, -c("N")] %>%
    melt(measure.vars = c("absolute_impact", "relative_impact", "NNV", "absolute_impact_p10000")) %>%
    .[, time := upper_time_limit] %>%
    return()
}) %>% rbindlist()
incidence_impact

incidence_impact_summarized = incidence_impact %>%
  .[, .(med = median(value),
      low50=quantile(value, 0.25), high50=quantile(value, 0.75),
      low95=quantile(value, 0.025), high95=quantile(value, 0.975)),
  by=c("age_group", "vaccine_strategy", "variable", "time", "routine_coverage")]

table_impact = incidence_impact_summarized %>%
  dcast(routine_coverage+age_group+vaccine_strategy+time~variable,
        value.var = c("med", "low95", "high95")) %>%
        .[, .(format = sprintf("%s (%s - %s, %s%%)",
                               round(med_absolute_impact_p10000),
                               round(low95_absolute_impact_p10000),
                               round(high95_absolute_impact_p10000),
                               round(med_relative_impact * 100))),
         by = c("age_group", "vaccine_strategy", "time", "routine_coverage")] %>%
  .[, time := factor(time, c(0.5, 1, 2, 3), c("0-6m", "0-1y", "0-2y", "0-3y"))] %>%
  .[, routine_coverage := factor(routine_coverage, names(routine_coverage_labels), routine_coverage_labels)] %>%
  dcast(age_group+time+vaccine_strategy~routine_coverage, value.var = "format") %>%
  merge(vaccination_strategies_plot_settings[, c("model_name", "name")], by.x = "vaccine_strategy", by.y = "model_name") %>%
  .[, vaccine_strategy := name] %>%
  .[, -"name"] %>%
  .[order(-age_group, time, vaccine_strategy)]
fwrite(table_impact, sprintf("%s/table_impact_routine.csv", OUTPUT_RESULTSFOLDER))

formatNumber = function(x){
  if(x > 1000) x = round(x/100) * 100
  format(round(x), big.mark = ",")
}

table_NNV = incidence_impact_summarized[variable == "NNV"] %>%
  .[, .(format = sprintf("%s (%s - %s)",
                         formatNumber(med),
                         formatNumber(low95),
                         formatNumber(high95))),
    by = c("age_group", "vaccine_strategy", "time", "routine_coverage")] %>%
  .[, time := factor(time, c(0.5, 1, 2, 3), c("0-6m", "0-1y", "0-2y", "0-3y"))] %>%
  .[, routine_coverage := factor(routine_coverage, names(routine_coverage_labels), routine_coverage_labels)] %>%
  dcast(age_group+time+vaccine_strategy~routine_coverage, value.var = "format") %>%
  merge(vaccination_strategies_plot_settings[, c("model_name", "name")], by.x = "vaccine_strategy", by.y = "model_name") %>%
  .[, vaccine_strategy := name] %>%
  .[, -"name"] %>%
  .[order(-age_group, time, vaccine_strategy)]
fwrite(table_NNV, sprintf("%s/table_NNV_routine.csv", OUTPUT_RESULTSFOLDER))